import re
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from general_utilities.association_resources import define_field_names_from_tarball_prefix, \
    bgzip_and_tabix
from general_utilities.import_utils.import_lib import download_bgen_file
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner


class SAIGERunner(ToolRunner):

    def run_tool(self) -> None:

        # Prep bgen files for a run:
        self._logger.info("Downloading and filtering raw bgen files")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A SAIGE bgen thread failed',
                                       incrementor=10,
                                       thread_factor=4)

        for chromosome in self._association_pack.bgen_dict:
            # This makes use of a utility class from AssociationResources since bgen filtering/processing is
            # IDENTICAL to that done for BOLT.
            thread_utility.launch_job(class_type=download_bgen_file,
                                      chrom_bgen_index=self._association_pack.bgen_dict[chromosome]
                                      )
        thread_utility.collect_futures()

        # 1. Run SAIGE step one without parallelisation
        self._logger.info("Running SAIGE step 1...")
        self._outputs.append(self._saige_step_one())

        # 2. Run SAIGE step two WITH parallelisation by chromosome
        self._logger.info("Running SAIGE step 2...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A SAIGE thread failed',
                                       incrementor=10,
                                       thread_factor=1)

        for chromosome in self._association_pack.bgen_dict:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                # changing this to search for BOLT files because we are now using these
                # as our .bgen inputs
                if Path(chromosome + ".bgen").exists():
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
        self._outputs.append(saige_step2_gene_log)

        # 4. Process final results
        self._logger.info("Processing final SAIGE output...")
        self._outputs.extend(self._annotate_saige_output(completed_gene_tables))

    # Run rare variant association testing using SAIGE-GENE
    def _saige_step_one(self) -> Path:

        # See the README.md for more information on these parameters
        # Just to note – I previously tried to implement the method that includes variance ratio estimation. However,
        # there are too few low MAC variants in the genotype files to perform this step accurately. The SAIGE
        # documentation includes this step, but I am very unsure how it works...
        cmd = f'step1_fitNULLGLMM.R ' \
              f'--phenoFile=/test/{self._association_pack.final_covariates} ' \
              f'--phenoCol={self._association_pack.pheno_names[0]} ' \
              f'--isCovariateTransform=FALSE ' \
              f'--sampleIDColinphenoFile=IID ' \
              f'--outputPrefix=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT ' \
              f'--sparseGRMFile=/test/{self._association_pack.sparse_grm.name} ' \
              f'--sparseGRMSampleIDFile=/test/{self._association_pack.sparse_grm_sample.name} ' \
              f'--nThreads={self._association_pack.threads} ' \
              f'--LOCO=FALSE ' \
              f'--skipModelFitting=FALSE ' \
              f'--useSparseGRMtoFitNULL=TRUE ' \
              f'--skipVarianceRatioEstimation=TRUE '
        if self._association_pack.is_binary:
            cmd = cmd + f'--traitType=binary '
        else:
            cmd = cmd + f'--traitType=quantitative '

        # Manage covariates
        if self._association_pack.ignore_base_covariates:
            all_covariates = []
            cat_covars = []
        else:
            all_covariates = [f'PC{PC}' for PC in range(1, 11)] + ['age', 'age_squared', 'wes_batch']
            if self._association_pack.sex == 2:
                all_covariates.append('sex')
            cat_covars = ['wes_batch']

        all_covariates.extend(self._association_pack.found_quantitative_covariates)
        all_covariates.extend(self._association_pack.found_categorical_covariates)
        cat_covars.extend(self._association_pack.found_categorical_covariates)

        if len(all_covariates) > 0:
            cmd = cmd + f'--covarColList=' + ','.join(all_covariates) + ' '
        if len(cat_covars) > 0:
            cmd = cmd + f'--qCovarColList=' + ','.join(cat_covars) + ' '

        saige_log_file = Path(f'{self._output_prefix}.SAIGE_step1.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=saige_log_file, print_cmd=True)

        return saige_log_file

    # This is a helper function to parallelise SAIGE step 2 by chromosome
    # This returns the tarball_prefix and chromosome number to make it easier to generate output
    def _saige_step_two(self, tarball_prefix: str, chromosome: str) -> Tuple[str, str, Path]:
        """
        Run SAIGE step 2 for a given chromosome.
        :param tarball_prefix: prefix for the tarball file (input)
        :param chromosome: chromosome / chunk to run SAIGE on (input)
        :return: tarball_prefix, chromosome, saige_log_file
        """

        # ACTION - we are going to run this function with the new Duat data
        # 1. run as it is now
        # 2. run with the addition of a filtered sample file and a flag (if exists) that can be used to filter
        # 3. run with a plink filtering command (worst case scenario) to filter the bgen file by sample inclusion

        # chromsomes should be stripped of the extras
        chromosome_num = re.match(r'chr(\d+)_', chromosome).group(1)

        # See the README.md for more information on these parameters
        cmd = f'step2_SPAtests.R ' \
              f'--bgenFile=/test/{chromosome}.bgen ' \
              f'--bgenFileIndex=/test/{chromosome}.bgen.bgi ' \
              f'--sampleFile=/test/{chromosome}.sample ' \
              f'--AlleleOrder=ref-first ' \
              f'--GMMATmodelFile=/test/{self._association_pack.pheno_names[0]}.SAIGE_OUT.rda ' \
              f'--sparseGRMFile=/test/{self._association_pack.sparse_grm.name} ' \
              f'--sparseGRMSampleIDFile=/test/{self._association_pack.sparse_grm_sample.name} ' \
              f'--LOCO=FALSE ' \
              f'--SAIGEOutputFile=/test/{tarball_prefix}.{chromosome}.SAIGE_OUT.SAIGE.gene.txt ' \
              f'--groupFile=/test/{tarball_prefix}.{chromosome}.SAIGE.groupFile.txt ' \
              f'--is_output_moreDetails=TRUE ' \
              f'--maxMAF_in_groupTest=0.5 ' \
              f'--maxMissing=1 ' \
              f'--chrom={chromosome_num} ' \
              f'--annotation_in_groupTest=foo '

        if self._association_pack.is_binary:
            cmd = cmd + '--is_Firth_beta=TRUE'

        saige_log_file = Path(f'{tarball_prefix}.{chromosome}.SAIGE_step2.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=saige_log_file)

        return tarball_prefix, chromosome, saige_log_file

    def _process_saige_output(self, tarball_prefix: str, chromosome: str) -> pd.DataFrame:
        # Load the raw table
        saige_table = pd.read_csv(tarball_prefix + "." + chromosome + ".SAIGE_OUT.SAIGE.gene.txt", sep='\t')
        saige_table = saige_table.rename(columns={'Region': 'ENST'})
        saige_table = saige_table.drop(columns=['Group', 'max_MAF'])

        # Get column names for Mask/MAF information if possible
        saige_table = define_field_names_from_tarball_prefix(tarball_prefix, saige_table)

        return saige_table

    def _annotate_saige_output(self, completed_gene_tables: list) -> List[Path]:
        # Create an output array
        plot_dir = Path(f'{self._output_prefix}_plots/')  # Path to store plots
        plot_dir.mkdir()
        outputs = [plot_dir]

        saige_table = pd.concat(completed_gene_tables)

        # Now merge the transcripts table into the gene table to add annotation and write
        saige_table = pd.merge(self._transcripts_table, saige_table, on='ENST', how="left")
        saige_path = Path(f'{self._output_prefix}.genes.SAIGE.stats.tsv')
        with saige_path.open('w') as gene_out:

            # Sort just in case
            saige_table = saige_table.sort_values(by=['chrom', 'start', 'end'])

            saige_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')

            # Make Manhattan plots
            for mask in saige_table['MASK'].value_counts().index:

                for maf in saige_table['MAF'].value_counts().index:
                    # To note on the below: I use SYMBOL for the id_column parameter below because ENST is the
                    # index and I don't currently have a way to pass the index through to the Plotter methods...
                    manhattan_plotter = ManhattanPlotter(self._association_pack.cmd_executor,
                                                         saige_table.query(f'MASK == "{mask}" & MAF == "{maf}"'),
                                                         chrom_column='chrom', pos_column='start',
                                                         alt_column=None,
                                                         id_column='ENST',
                                                         p_column='Pvalue',
                                                         csq_column='MASK',
                                                         maf_column='MAC', gene_symbol_column='SYMBOL',
                                                         clumping_distance=1,
                                                         maf_cutoff=30,
                                                         sig_threshold=1E-6)

                    manhattan_plotter.plot()[0].rename(plot_dir / f'{mask}.{maf}.genes.SAIGE.png')

        # And bgzip and tabix...
        outputs.extend(bgzip_and_tabix(saige_path, skip_row=1, sequence_row=2, begin_row=3, end_row=4))

        return outputs
