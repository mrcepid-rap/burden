import pandas as pd

from typing import List, Tuple
from pathlib import Path

from burden.tool_runners.tool_runner import ToolRunner
from general_utilities.association_resources import get_chromosomes, process_bgen_file, \
    define_field_names_from_tarball_prefix, build_transcript_table, bgzip_and_tabix
from general_utilities.job_management.thread_utility import ThreadUtility


class SAIGERunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Run SAIGE step one without parallelisation
        self._logger.info("Running SAIGE step 1...")
        self._outputs.append(self._saige_step_one())

        # 2. Run SAIGE step two WITH parallelisation by chromosome
        self._logger.info("Running SAIGE step 2...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A SAIGE thread failed', 
                                       incrementor=10, 
                                       thread_factor=1)
        
        for chromosome in get_chromosomes():
            for tarball_prefix in self._association_pack.tarball_prefixes:
                if Path(tarball_prefix + "." + chromosome + ".SAIGE.bcf").exists():
                    self._prep_group_file(tarball_prefix, chromosome)
                    thread_utility.launch_job(class_type=self._saige_step_two,
                                              tarball_prefix=tarball_prefix,
                                              chromosome=chromosome)

        # 3. Gather preliminary results
        self._logger.info("Gathering SAIGE mask-based results...")
        completed_gene_tables = []
        saige_step2_gene_log = Path(f'{self._output_prefix}.SAIGE_step2.log')
        with saige_step2_gene_log.open('w') as saige_step2_genes_writer:
            for result in thread_utility:
                tarball_prefix, finished_chromosome, current_log = result
                completed_gene_tables.append(self._process_saige_output(tarball_prefix, finished_chromosome))

                # Write a header for each file
                saige_step2_genes_writer.write(f'{tarball_prefix + "-" + finished_chromosome:{"-"}^{50}}')
                with current_log.open('r') as current_log_reader:
                    for line in current_log_reader:
                        saige_step2_genes_writer.write(line)
        completed_gene_tables.append(saige_step2_gene_log)

        # 4. Run per-marker tests, if requested
        completed_marker_chromosomes = []
        if self._association_pack.run_marker_tests:
            self._logger.info("Running per-marker tests...")
            thread_utility = ThreadUtility(self._association_pack.threads,
                                           error_message='A SAIGE marker thread failed',
                                           incrementor=1,
                                           thread_factor=4)

            for chromosome in get_chromosomes():
                thread_utility.launch_job(class_type=self._saige_marker_run,
                                          chromosome=chromosome)
                completed_marker_chromosomes.append(chromosome)

            saige_step2_markers_log = Path(f'{self._output_prefix}.SAIGE_step2.log')
            with saige_step2_markers_log.open('w') as saige_step2_markers_writer:
                for result in thread_utility:
                    finished_chromosome, current_log = result
                    saige_step2_markers_writer.write(f'{finished_chromosome + ".log":{"-"}^{50}}')
                    with current_log.open('r') as current_log_reader:
                        for line in current_log_reader:
                            saige_step2_markers_writer.write(line)
            self._outputs.append(saige_step2_markers_log)

        # 5. Process final results
        self._logger.info("Processing final SAIGE output...")
        self._outputs.extend(self._annotate_saige_output(completed_gene_tables, completed_marker_chromosomes))

    # Run rare variant association testing using SAIGE-GENE
    def _saige_step_one(self) -> Path:

        # See the README.md for more information on these parameters
        # Just to note – I previously tried to implement the method that includes variance ratio estimation. However,
        # there are too few low MAC variants in the genotype files to perform this step accurately. The SAIGE
        # documentation includes this step, but I am very unsure how it works...
        cmd = 'step1_fitNULLGLMM.R ' \
                    '--phenoFile=/test/phenotypes_covariates.formatted.txt ' \
                    f'--phenoCol={self._association_pack.pheno_names[0]} ' \
                    '--isCovariateTransform=FALSE ' \
                    '--sampleIDColinphenoFile=IID ' \
                    f'--outputPrefix=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT ' \
                    '--sparseGRMFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx ' \
                    '--sparseGRMSampleIDFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt ' \
                    f'--nThreads={str(self._association_pack.threads)} ' \
                    '--LOCO=FALSE ' \
                    '--skipModelFitting=FALSE ' \
                    '--useSparseGRMtoFitNULL=TRUE ' \
                    '--skipVarianceRatioEstimation=TRUE '
        if self._association_pack.is_binary:
            cmd = cmd + '--traitType=binary '
        else:
            cmd = cmd + '--traitType=quantitative '

        # Manage covariates
        if self._association_pack.sex == 2:
            default_covars = ['PC' + str(x) for x in range(1, 11)] + ['age', 'age_squared', 'sex', 'wes_batch']
        else:
            default_covars = ['PC' + str(x) for x in range(1, 11)] + ['age', 'age_squared', 'wes_batch']
        all_covariates = [','.join(default_covars)]  # A list to manage appending additional covariates
        if len(self._association_pack.found_quantitative_covariates) > 0:
            all_covariates.append(','.join(self._association_pack.found_quantitative_covariates))
        if len(self._association_pack.found_categorical_covariates) > 0:
            cat_covars_join = ','.join(self._association_pack.found_categorical_covariates)
            all_covariates.append(cat_covars_join)
            cmd = cmd + '--qCovarColList=wes_batch,' + cat_covars_join + ' '
        else:
            cmd = cmd + '--qCovarColList=wes_batch '
        cmd = cmd + '--covarColList=' + ','.join(all_covariates)

        saige_log_file = Path(f'{self._output_prefix}.SAIGE_step1.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=saige_log_file, print_cmd=True)

        return saige_log_file

    # This exists for a very stupid reason – they _heavily_ modified the groupFile for v1.0 and I haven't gone back
    # to change how this file is made in 'collapse variants'
    @staticmethod
    def _prep_group_file(tarball_prefix: str, chromosome: str) -> None:

        with open(tarball_prefix + '.' + chromosome + '.SAIGE.groupFile.txt', 'r') as group_file:
            modified_group = open(tarball_prefix + '.' + chromosome + '.SAIGE_v1.0.groupFile.txt', 'w')
            for line in group_file:
                data = line.rstrip().split('\t')
                mod_data = [data[0], 'var']
                found = set()
                for var in data[1:]:
                    var = var.translate(str.maketrans('_/', '::'))
                    if var not in found:
                        mod_data.append(var)
                        found.add(var)
                modified_group.writelines(' '.join(mod_data) + "\n")
                mod_annote = [data[0], 'anno']
                for i in range(2, len(mod_data)):
                    mod_annote.append('foo')
                modified_group.write(' '.join(mod_annote) + "\n")
            group_file.close()
            modified_group.close()

    # This is a helper function to parallelise SAIGE step 2 by chromosome
    # This returns the tarball_prefix and chromosome number to make it easier to generate output
    def _saige_step_two(self, tarball_prefix: str, chromosome: str) -> Tuple[str, str, Path]:

        cmd = f'bcftools view --threads 1 -S /test/SAMPLES_Include.txt -Ob ' \
              f'-o /test/{tarball_prefix}.{chromosome}.saige_input.bcf ' \
              f'/test/{tarball_prefix}.{chromosome}.SAIGE.bcf'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        cmd = f'bcftools index --threads 1 /test/{tarball_prefix}.{chromosome}.saige_input.bcf'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)
        
        # See the README.md for more information on these parameters
        cmd = f'step2_SPAtests.R ' \
              f'--vcfFile=/test/{tarball_prefix}.{chromosome}.saige_input.bcf ' \
              f'--vcfField=GT ' \
              f'--GMMATmodelFile=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT.rda ' \
              f'--sparseGRMFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx ' \
              f'--sparseGRMSampleIDFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt ' \
              f'--LOCO=FALSE ' \
              f'--SAIGEOutputFile=/test/{tarball_prefix}.{chromosome}.SAIGE_OUT.SAIGE.gene.txt ' \
              f'--groupFile=/test/{tarball_prefix}.{chromosome}.SAIGE_v1.0.groupFile.txt ' \
              f'--is_output_moreDetails=TRUE ' \
              f'--maxMAF_in_groupTest=0.5 ' \
              f'--maxMissing=1 ' \
              f'--chrom={chromosome} ' \
              f'--annotation_in_groupTest=foo '

        if self._association_pack.is_binary:
            cmd = cmd + '--is_Firth_beta=TRUE'

        saige_log_file = Path(f'{tarball_prefix}.{chromosome}.SAIGE_step2.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=saige_log_file)

        return tarball_prefix, chromosome, saige_log_file

    def _saige_marker_run(self, chromosome: str) -> Tuple[str, Path]:

        process_bgen_file(self._association_pack.bgen_dict[chromosome], chromosome)

        cmd = 'step2_SPAtests.R ' \
              f'--bgenFile=/test/{chromosome}.markers.bgen ' \
              f'--bgenFileIndex=/test/{chromosome}.markers.bgen.bgi ' \
              f'--sampleFile=/test/{chromosome}.markers.sample ' \
              f'--GMMATmodelFile=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT.rda ' \
              '--sparseGRMFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx ' \
              '--sparseGRMSampleIDFile=/test/genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt ' \
              f'--SAIGEOutputFile=/test/{chromosome}.SAIGE_OUT.SAIGE.markers.txt ' \
              '--LOCO=FALSE ' \
              '--is_output_moreDetails=TRUE ' \
              '--maxMissing=1 '
        if self._association_pack.is_binary:
            cmd = cmd + '--is_Firth_beta=TRUE'

        saige_log_file = Path(f'{chromosome}.SAIGE_markers.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=saige_log_file)

        return chromosome, saige_log_file

    def _process_saige_output(self, tarball_prefix: str, chromosome: str) -> pd.DataFrame:

        # Load the raw table
        saige_table = pd.read_csv(tarball_prefix + "." + chromosome + ".SAIGE_OUT.SAIGE.gene.txt", sep='\t')
        saige_table = saige_table.rename(columns={'Region': 'ENST'})
        saige_table = saige_table.drop(columns=['Group', 'max_MAF'])

        # Get column names for Mask/MAF information if possible
        saige_table = define_field_names_from_tarball_prefix(tarball_prefix, saige_table)

        return saige_table

    def _annotate_saige_output(self, completed_gene_tables: list, completed_marker_chromosomes: list) -> List[Path]:

        # Create an output array
        outputs = []

        saige_table = pd.concat(completed_gene_tables)

        # Now process the gene table into a useable format:
        # First read in the transcripts file
        transcripts_table = build_transcript_table()

        # Now merge the transcripts table into the gene table to add annotation and write
        saige_table = pd.merge(transcripts_table, saige_table, on='ENST', how="left")
        saige_path = Path(f'{self._output_prefix}.genes.SAIGE.stats.tsv')
        with saige_path.open('w') as gene_out:

            # Sort just in case
            saige_table = saige_table.sort_values(by=['chrom', 'start', 'end'])

            saige_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')

        # And bgzip and tabix...
        outputs.extend(bgzip_and_tabix(saige_path, skip_row=1, sequence_row=2, begin_row=3, end_row=4))

        if self._association_pack.run_marker_tests:

            variant_index = []
            saige_table_marker = []
            # Open all chromosome indicies and load them into a list and append them together
            for chromosome in completed_marker_chromosomes:
                variant_index.append(pd.read_csv(f'filtered_bgen/{chromosome}.filtered.vep.tsv.gz',
                                                 sep="\t",
                                                 dtype={'SIFT': str, 'POLYPHEN': str}))
                saige_table_marker.append(pd.read_csv(chromosome + ".SAIGE_OUT.SAIGE.markers.txt", sep="\t"))

            variant_index = pd.concat(variant_index)
            variant_index = variant_index.set_index('varID')

            saige_table_marker = pd.concat(saige_table_marker)

            # For markers, we can use the SNP ID column to get what we need
            saige_table_marker = saige_table_marker.rename(columns={'MarkerID': 'varID',
                                                                    'AC_Allele2': 'SAIGE_AC',
                                                                    'AF_Allele2': 'SAIGE_MAF'})
            saige_table_marker = saige_table_marker.drop(columns=['CHR', 'POS', 'Allele1', 'Allele2', 'MissingRate'])

            saige_table_marker = pd.merge(variant_index, saige_table_marker, on='varID', how="left")
            saige_marker_path = Path(f'{self._output_prefix}.markers.SAIGE.stats.tsv')
            with saige_marker_path.open('w') as marker_out:
                # Sort by chrom/pos just to be sure...
                saige_table_marker = saige_table_marker.sort_values(by=['CHROM', 'POS'])

                saige_table_marker.to_csv(path_or_buf=marker_out, index=False, sep="\t", na_rep='NA')

            # And bgzip and tabix...
            outputs.extend(bgzip_and_tabix(saige_marker_path, skip_row=1, sequence_row=2, begin_row=3))

        return outputs
