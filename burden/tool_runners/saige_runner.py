import re
from pathlib import Path
from typing import List, Tuple, Dict, Any

import dxpy
import pandas as pd
from general_utilities.association_resources import define_field_names_from_tarball_prefix, \
    bgzip_and_tabix, LOGGER
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner


class SAIGERunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Run SAIGE step one without parallelisation
        self._logger.info("Running SAIGE step 1...")
        self._outputs.append(self._saige_step_one())

        # 2. Run SAIGE step two WITH parallelisation by chromosome
        self._logger.info("Running SAIGE step 2 with subjobs...")

        all_step2_outputs = self._multithread_step2()

        # 3. Gather preliminary results
        self._logger.info("Gathering SAIGE mask-based results...")
        completed_gene_tables = []
        saige_step2_gene_log = Path(f'{self._output_prefix}.SAIGE_step2.log')
        with saige_step2_gene_log.open('w') as saige_step2_genes_writer:
            for result in all_step2_outputs:
                tarball_prefix = result['tarball_prefix']
                finished_chromosome = result['finished_chromosome']
                current_log = InputFileHandler(result['saige_log_file']).get_file_handle()

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

        # SAIGE will complain if the BGEN contains samples that are not in the phenotype file, so let's make sure
        # we subset just in case
        first_sample_key = next(iter(self._association_pack.bgen_dict))
        sample = pd.read_csv(
            self._association_pack.bgen_dict[first_sample_key]['sample'].get_file_handle(),
            delim_whitespace=True,
            dtype=str
        )
        print(sample)
        phenotype = pd.read_csv(self._association_pack.final_covariates, sep=' ', dtype=str)
        print(phenotype)
        # Subset samples where the ID is in the phenotype file
        subset = phenotype[phenotype.iloc[:, 0].isin(sample.iloc[:, 0])]
        # the second row must be 0 / 0 / D
        print(subset)
        # Save the result
        phenofile = Path("phenotype_subset_sample.txt")
        subset.to_csv(phenofile, sep='\t', index=False)

        # See the README.md for more information on these parameters
        # Just to note – I previously tried to implement the method that includes variance ratio estimation. However,
        # there are too few low MAC variants in the genotype files to perform this step accurately. The SAIGE
        # documentation includes this step, but I am very unsure how it works...
        cmd = f'step1_fitNULLGLMM.R ' \
              f'--phenoFile={phenofile} ' \
              f'--phenoCol={self._association_pack.pheno_names[0]} ' \
              f'--isCovariateTransform=FALSE ' \
              f'--sampleIDColinphenoFile=IID ' \
              f'--outputPrefix={self._association_pack.pheno_names[0]}.SAIGE_OUT ' \
              f'--sparseGRMFile={self._association_pack.sparse_grm} ' \
              f'--sparseGRMSampleIDFile={self._association_pack.sparse_grm_sample} ' \
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
            all_covariates = [f'PC{PC}' for PC in range(1, 11)] + ['age', 'age_squared', 'batch']
            if self._association_pack.sex == 2:
                all_covariates.append('sex')
            cat_covars = ['batch']

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

    def _multithread_step2(self) -> List[Dict[str, Any]]:
        """
        A wrapper function to allow for multithreading of SAIGE step 2 by chromosome.

        :return: A list of dictionaries containing the outputs from each subjob
        """

        # set the launcher
        launcher = joblauncher_factory()

        # set the exporter
        exporter = ExportFileHandler(delete_on_upload=False)

        for chromosome in self._association_pack.bgen_dict:

            # make a list of the setlist files for this chromosome
            group_files = list(Path('.').glob(f'*.{chromosome}.SAIGE.groupFile.txt'))

            # export the files for each subjob
            gmmatmodelfile = exporter.export_files(f"{self._association_pack.pheno_names[0]}.SAIGE_OUT.rda")
            sparsegrmfile = exporter.export_files(f"{self._association_pack.sparse_grm}")
            sparsegrmsampleidfile = exporter.export_files(f"{self._association_pack.sparse_grm_sample}")

            group_files = [exporter.export_files(gf) for gf in group_files]

            launcher.launch_job(
                function=run_saige_step_two,
                inputs={
                    'bgen_file': self._association_pack.bgen_dict[chromosome]['bgen'].get_input_str(),
                    'bgen_index': self._association_pack.bgen_dict[chromosome]['index'].get_input_str(),
                    'sample_file':self._association_pack.bgen_dict[chromosome]['sample'].get_input_str(),
                    'chromosome': chromosome,
                    "tarball_prefixes": self._association_pack.tarball_prefixes,
                    'gmmatmodelfile': gmmatmodelfile,
                    'sparsegrmfile': sparsegrmfile,
                    'sparsegrmsampleidfile': sparsegrmsampleidfile,
                    'group_files': group_files,
                    'is_binary': self._association_pack.is_binary
                },
                outputs=['output'],
            )
        launcher.submit_and_monitor()

        step2_outputs = []
        for result in launcher:
            # result["output"] is already a list of dicts
            for r in result["output"]:
                # download subjob outputs to local machine
                InputFileHandler(r["saige_log_file"], download_now=True)
                InputFileHandler(r["saige_output"], download_now=True)
                step2_outputs.append(r)

        return step2_outputs

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


@dxpy.entry_point('run_saige_step_two')
def run_saige_step_two(bgen_file: str, bgen_index: str, sample_file: str,
                       chromosome: str, tarball_prefixes: List[str], gmmatmodelfile: str,
                       sparsegrmfile: str, sparsegrmsampleidfile: str, group_files: List[str],
                       is_binary: bool) -> Dict[str, Any]:
    """
    A wrapper function to allow for multithreading of SAIGE step 2 by chromosome.

    :param bgen_file: The bgen file to use
    :param bgen_index: The bgen index file to use
    :param sample_file: The sample file to use
    :param chromosome: The chromosome / chunk to run SAIGE on
    :param tarball_prefixes: The tarball prefixes to run SAIGE on
    :param gmmatmodelfile: The GMMAT model file from step 1
    :param sparsegrmfile: The sparse GRM file to use
    :param sparsegrmsampleidfile: The sparse GRM sample ID file to use
    :param group_files: The group files to use
    :param is_binary: Is the phenotype binary?

    :return:
    """
    # we are now working per chunk
    # Get all the files we need
    bgen_file = InputFileHandler(bgen_file).get_file_handle()
    bgen_index = InputFileHandler(bgen_index).get_file_handle()
    sample_file = InputFileHandler(sample_file).get_file_handle()
    gmmatmodelfile = InputFileHandler(gmmatmodelfile).get_file_handle()
    sparsegrmfile = InputFileHandler(sparsegrmfile).get_file_handle()
    sparsegrmsampleidfile = InputFileHandler(sparsegrmsampleidfile).get_file_handle()
    for group_file in group_files:
        InputFileHandler(group_file, download_now=True)

    # 4. Run step 2 of SAIGE
    LOGGER.info("Running SAIGE step 2")
    thread_utility = ThreadUtility(thread_factor=1)

    for tarball_prefix in tarball_prefixes:

        LOGGER.info(f"Running for the mask {tarball_prefix}")

        # if Path(f'{tarball_prefix}.{chromosome}.REGENIE.annotationFile.tsv').exists():
        if Path(f'{tarball_prefix}.{chromosome}.SAIGE.groupFile.txt').exists():
            thread_utility.launch_job(function=saige_step_two,
                                      inputs={
                                          "tarball_prefix": tarball_prefix,
                                          "chromosome": chromosome,
                                          "bgen_file": bgen_file,
                                          "bgen_index": bgen_index,
                                          "sample_file": sample_file,
                                          "gmmatmodelfile": gmmatmodelfile,
                                          "sparsegrmfile": sparsegrmfile,
                                          "sparsegrmsampleidfile": sparsegrmsampleidfile,
                                          "is_binary": is_binary
                                      },
                                      outputs=[
                                          "tarball_prefix",
                                          "chromosome",
                                          "saige_log_file",
                                          "saige_output"
                                      ]
                                      )
    thread_utility.submit_and_monitor()

    # collect results from thread_utility and export files
    exporter = ExportFileHandler()
    output = []
    for result in thread_utility:
        output.append({
            "tarball_prefix": result["tarball_prefix"],
            "finished_chromosome": result["chromosome"],
            "saige_log_file": exporter.export_files(result["saige_log_file"]),
            "saige_output": exporter.export_files(result["saige_output"])
        })

    return {
        "output": output
    }


# This is a helper function to parallelise SAIGE step 2 by chromosome
# This returns the tarball_prefix and chromosome number to make it easier to generate output
def saige_step_two(tarball_prefix: str, chromosome: str, bgen_file, bgen_index, sample_file, gmmatmodelfile,
                   sparsegrmfile, sparsegrmsampleidfile, is_binary) -> Tuple[str, str, Path, Path]:
    """
    Run SAIGE step 2 for a given chromosome.
    :param tarball_prefix: prefix for the tarball file (input)
    :param chromosome: chromosome / chunk to run SAIGE on (input)
    :param bgen_file: The bgen file to use
    :param bgen_index: The bgen index file to use
    :param sample_file: The sample file to use
    :param gmmatmodelfile: The GMMAT model file from step 1
    :param sparsegrmfile: The sparse GRM file to use
    :param sparsegrmsampleidfile: The sparse GRM sample ID file to use
    :param is_binary: Is the phenotype binary?

    :return: tarball_prefix, chromosome, saige_log_file
    """

    cmd_exec = build_default_command_executor()

    # ACTION - we are going to run this function with the new Duat data
    # 1. run as it is now
    # 2. run with the addition of a filtered sample file and a flag (if exists) that can be used to filter
    # 3. run with a plink filtering command (worst case scenario) to filter the bgen file by sample inclusion

    # chromsomes should be stripped of
    if chromosome.startswith("chr"):
        chromosome_num = re.match(r'chr(\d+)_', chromosome).group(1)
    else:
        chromosome_num = chromosome

    sample_df = pd.read_csv(sample_file, delim_whitespace=True, low_memory=False)
    print(sample_df.head())
    # if the columns are ID | missing | sex
    if list(sample_df.columns) == ['ID', 'missing', 'sex']:
        LOGGER.info(f"Formatting sample file for REGENIE (likely from WES data)")
        # Duplicate ID into ID_1 and ID_2
        sample_df['ID_1'] = sample_df['ID']
        sample_df['ID_2'] = sample_df['ID']
        # Reorder columns
        sample_df = sample_df[['ID_1', 'ID_2', 'missing', 'sex']]
        print(sample_df)
        sample_df.to_csv(Path(sample_file), sep="\t", index=False)

    # See the README.md for more information on these parameters
    cmd = f'step2_SPAtests.R ' \
          f'--bgenFile={bgen_file} ' \
          f'--bgenFileIndex={bgen_index} ' \
          f'--sampleFile={sample_file} ' \
          f'--AlleleOrder=ref-first ' \
          f'--GMMATmodelFile={gmmatmodelfile} ' \
          f'--sparseGRMFile={sparsegrmfile} ' \
          f'--sparseGRMSampleIDFile={sparsegrmsampleidfile} ' \
          f'--LOCO=FALSE ' \
          f'--SAIGEOutputFile={tarball_prefix}.{chromosome}.SAIGE_OUT.SAIGE.gene.txt ' \
          f'--groupFile={tarball_prefix}.{chromosome}.SAIGE.groupFile.txt ' \
          f'--is_output_moreDetails=TRUE ' \
          f'--maxMAF_in_groupTest=0.5 ' \
          f'--maxMissing=1 ' \
          f'--chrom={chromosome_num} ' \
          f'--annotation_in_groupTest=foo '

    if is_binary:
        cmd = cmd + '--is_Firth_beta=TRUE'

    saige_log_file = Path(f'*.SAIGE_step2.log')
    cmd_exec.run_cmd_on_docker(cmd, stdout_file=saige_log_file, print_cmd=True)

    saige_output = Path(f'{tarball_prefix}.{chromosome}.SAIGE_OUT.SAIGE.gene.txt')

    return tarball_prefix, chromosome, saige_log_file, saige_output
