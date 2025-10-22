import re
from pathlib import Path
from typing import List, Any, Tuple, Dict

import dxpy
import pandas as pd
from general_utilities.association_resources import define_covariate_string, \
    bgzip_and_tabix, define_field_names_from_tarball_prefix
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import LOGGER
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner


# TODO: Implement multi-phenotype testing for REGENIE
class REGENIERunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Run step 1 of regenie
        # if step one has already been run (if the files exist and are not empty)
        if not self._association_pack.has_regenie_step_one:
            self._logger.info("Running REGENIE step 1")
            regenie_step1_log = self._run_regenie_step_one()
            # Add the step1 files to output, so we can use later if need-be:
            self._outputs.extend([Path('fit_out_pred.list'),
                                  Path('fit_out_1.loco'),
                                  regenie_step1_log])
        else:
            self._logger.info("REGENIE step 1 files already exist, skipping...")

        # # 2. Prep bgen files for a run:
        self._logger.info("Preparing to launch subjobs for running REGENIE step 2 on separate VMs...")

        all_step2_outputs = self._multithread_step2()

        # Gather preliminary results from step 2:
        self._logger.info("Gathering REGENIE mask-based results...")
        completed_gene_tables = []
        regenie_step2_genes_log = Path(f'{self._output_prefix}.REGENIE_step2_genes.log')
        with regenie_step2_genes_log.open('w') as regenie_step2_genes_writer:
            for result in all_step2_outputs:
                tarball_prefix = result['tarball_prefix']
                finished_chromosome = result['finished_chromosome']
                current_log = InputFileHandler(result['current_log']).get_file_handle()

                completed_gene_tables.append(self._process_regenie_output(tarball_prefix,
                                                                          finished_chromosome))

                # Write a header for each logfile
                regenie_step2_genes_writer.write(f'{tarball_prefix + "-" + finished_chromosome:{"-"}^{50}}')

                with current_log.open('r') as current_log_reader:
                    for line in current_log_reader:
                        regenie_step2_genes_writer.write(line)

            self._outputs.append(regenie_step2_genes_log)

            # Always add predictions...
            self._outputs.append(Path('fit_out_pred.list'))
            self._outputs.extend(Path('./').glob('fit_out_*.loco'))

            # 5. Process outputs
            self._logger.info("Processing REGENIE outputs...")
            self._outputs.extend(self._annotate_regenie_output(completed_gene_tables=completed_gene_tables))

    def _run_regenie_step_one(self) -> Path:

        # And generate a SNP list for the --extract parameter of REGENIE, while considering SNPs from
        # the regenie_smaller_snps input parameter (if provided). plink2 order of operations:
        # 1. Select variants from --extract (if present)
        # 2. THEN filter based on max/min AC (mac/max-mac)
        max_mac = (self._sample_count * 2) - 100
        cmd = f"plink2 --bed {self._association_pack.genetic_filename}.bed " \
              f"--bim {self._association_pack.genetic_filename}.bim " \
              f"--fam {self._association_pack.genetic_filename}.fam " \
              f"--min-ac 100 " \
              f'--max-ac {str(max_mac)}' \
              f' --write-snplist ' \
              f'--out REGENIE_extract'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        with open('REGENIE_extract.snplist', 'r') as plink_out:
            for line in plink_out:
                found_snp_count = re.search('(\\d+) variants remaining after main filters', line)
                if found_snp_count is not None:
                    self._logger.info(f'Number of SNPs for REGENIE Step 1: {found_snp_count.group(1)}')
            plink_out.close()

        cmd = 'regenie ' \
              '--step 1 ' \
              f'--bed {self._association_pack.genetic_filename} ' \
              '--extract REGENIE_extract.snplist ' \
              f'--covarFile {self._association_pack.final_covariates} ' \
              f'--phenoFile {self._association_pack.final_covariates} ' \
              '--maxCatLevels 110 ' \
              '--bsize 1000 ' \
              '--out fit_out ' \
              f'--threads {str(self._association_pack.threads)} ' \
              f'--phenoCol {self._association_pack.pheno_names[0]} '

        cmd += define_covariate_string(self._association_pack.found_quantitative_covariates,
                                       self._association_pack.found_categorical_covariates,
                                       self._association_pack.is_binary,
                                       add_array=False,
                                       ignore_base=self._association_pack.ignore_base_covariates)

        regenie_log = Path(f'{self._output_prefix}.REGENIE_step1.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=regenie_log)

        return regenie_log

    def _multithread_step2(self) -> List[Dict[str, Any]]:
        """
        A function to run REGENIE step 2 on a single chromosome in a multithreaded manner

        :return: A list of dictionaries containing the tarball prefix, finished chromosome, and log file
        """

        # set the launcher
        launcher = joblauncher_factory()

        # files to include
        samples_include = Path('SAMPLES_Include.txt')
        fit_out_pred = Path('fit_out_pred.list')
        fit_out_loco = Path('fit_out_1.loco')

        # set the exporter
        exporter = ExportFileHandler(delete_on_upload=False)

        # print all files in the current directory
        for path in Path('.').glob('*'):
            LOGGER.info(f'Found file: {path}')

        print("chromosomes to run:", self._association_pack.bgen_dict.keys())
        print("masks to run:", self._association_pack.tarball_prefixes)

        for chromosome in self._association_pack.bgen_dict:

            # make a list of all annotation files for this chromosome
            anno_files = list(Path('.').glob(f'*.{chromosome}.REGENIE.annotationFile.txt'))
            if not anno_files:
                continue

            # make a list of the mask files for this chromosome
            mask_files = list(Path('.').glob(f'*.{chromosome}.REGENIE.maskfile.txt'))
            if not mask_files:
                continue

            # make a list of the setlist files for this chromosome
            setlist_files = list(Path('.').glob(f'*.{chromosome}.REGENIE.setListFile.txt'))
            if not setlist_files:
                continue

            # export the files to DX for each subjob
            samples_include = exporter.export_files(samples_include)
            fit_out_pred = exporter.export_files(fit_out_pred)
            fit_out_loco = exporter.export_files(fit_out_loco)
            phenotype_file = exporter.export_files(self._association_pack.final_covariates)
            anno_files = [exporter.export_files(af) for af in anno_files]
            mask_files = [exporter.export_files(mf) for mf in mask_files]
            setlist_files = [exporter.export_files(sf) for sf in setlist_files]

            launcher.launch_job(function=run_regenie_step2,
                                inputs={
                                    "bgen_file": self._association_pack.bgen_dict[chromosome]['bgen'].get_input_str(),
                                    "bgen_sample": self._association_pack.bgen_dict[chromosome][
                                        'sample'].get_input_str(),
                                    "bgen_index": self._association_pack.bgen_dict[chromosome]['index'].get_input_str(),
                                    "chromosome": chromosome,
                                    "tarball_prefixes": self._association_pack.tarball_prefixes,
                                    "samples_include": samples_include,
                                    "covariate_file": phenotype_file,
                                    "pheno_file": phenotype_file,
                                    "pheno_column": self._association_pack.pheno_names[0],
                                    "fit_out_pred": fit_out_pred,
                                    "fit_out_loco": fit_out_loco,
                                    "annotation_file": anno_files,
                                    "mask_file": mask_files,
                                    "setlist_file": setlist_files,
                                    "found_quantitative_covariates": self._association_pack.found_quantitative_covariates,
                                    "found_categorical_covariates": self._association_pack.found_categorical_covariates,
                                    "is_binary": self._association_pack.is_binary,
                                    "ignore_base_covariates": self._association_pack.ignore_base_covariates,
                                },
                                outputs=[
                                    "output"
                                ]
                                )
        launcher.submit_and_monitor()

        step2_outputs = []
        for result in launcher:
            # result["output"] is already a list of dicts
            for r in result["output"]:
                # download subjob outputs to local machine
                InputFileHandler(r["current_log"], download_now=True)
                InputFileHandler(r["regenie_output"], download_now=True)
                step2_outputs.append(r)

        return step2_outputs

    def _process_regenie_output(self, tarball_prefix: str, chromosome: str) -> pd.DataFrame:

        regenie_table = pd.read_csv(f'{tarball_prefix}.{chromosome}_{self._association_pack.pheno_names[0]}.regenie',
                                    sep=' ',
                                    comment='#')

        # And then should be able to split into 3 columns:
        regenie_table[['ENST', 'MASK', 'SUBSET']] = regenie_table['ID'].str.split('.', expand=True)

        # And only take the 'all' subset. Singleton is meaningless here
        regenie_table = regenie_table[regenie_table['SUBSET'] == 'all']

        # Convert -log10(p) back into a p. value:
        regenie_table['PVALUE'] = 10 ** (-1 * regenie_table['LOG10P'])

        # And finally drop columns we won't care about:
        regenie_table = regenie_table.drop(columns=['ID', 'CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'EXTRA',
                                                    'MASK', 'SUBSET', 'LOG10P'])

        # Get column names for Mask/MAF information if possible from the tarball name
        regenie_table = define_field_names_from_tarball_prefix(tarball_prefix, regenie_table)

        # Here we pull out all the various tests that were performed and make a separate table
        p_val_table = regenie_table.pivot(index='ENST', columns='TEST', values='PVALUE')
        p_val_table = p_val_table.drop(columns=['ADD'])

        # And then merge them back into the original table with just the 'additive' p. value
        regenie_table = regenie_table.set_index('ENST')
        regenie_table = regenie_table[regenie_table['TEST'] == 'ADD']
        regenie_table = regenie_table.merge(p_val_table, left_index=True, right_index=True)

        # And finally add an AC column
        regenie_table['AC'] = regenie_table['A1FREQ'] * (regenie_table['N'] * 2)
        regenie_table['AC'] = regenie_table['AC'].round().astype(int)

        return regenie_table

    def _annotate_regenie_output(self, completed_gene_tables: List[pd.DataFrame]) -> List[Path]:

        # Declare a list for storing outputs
        plot_dir = Path(f'{self._output_prefix}_plots/')  # Path to store plots
        plot_dir.mkdir(exist_ok=True)
        outputs = [plot_dir]
        regenie_table = pd.concat(completed_gene_tables)

        # Now merge the transcripts table into the gene table to add annotation and then write
        regenie_table = pd.merge(self._transcripts_table, regenie_table, left_index=True, right_index=True,
                                 how="left")

        regenie_gene_out = Path(f'{self._output_prefix}.genes.REGENIE.stats.tsv')
        with regenie_gene_out.open('w') as gene_out:

            # Reset the index and make sure chrom/start/end are first (for indexing)
            regenie_table.reset_index(inplace=True)
            columns = regenie_table.columns.tolist()

            # ENST should ALWAYS be in position 0, but move it to position 4 and slice the array so we don't have two
            # copies:

            columns = columns[1:]
            regenie_table = regenie_table[columns]

            # Sort just in case
            regenie_table = regenie_table.sort_values(by=['chrom', 'start', 'end'])

            for mask in regenie_table['MASK'].value_counts().index:

                for maf in regenie_table['MAF'].value_counts().index:
                    # To note on the below: I use SYMBOL for the id_column parameter below because ENST is the
                    # index and I don't currently have a way to pass the index through to the Plotter methods...
                    manhattan_plotter = ManhattanPlotter(self._association_pack.cmd_executor,
                                                         regenie_table.query(f'MASK == "{mask}" & MAF == "{maf}"'),
                                                         chrom_column='chrom', pos_column='start',
                                                         alt_column=None,
                                                         id_column='SYMBOL', p_column='PVALUE',
                                                         csq_column='MASK',
                                                         maf_column='A1FREQ', gene_symbol_column='SYMBOL',
                                                         clumping_distance=1,
                                                         maf_cutoff=30 / (regenie_table['N'].max() * 2),
                                                         sig_threshold=1E-6)

                    manhattan_plotter.plot()[0].rename(plot_dir / f'{mask}.{maf}.genes.REGENIE.png')
            # Write to disk
            regenie_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')

        # And bgzip and tabix...
        outputs.extend(
            bgzip_and_tabix(regenie_gene_out, comment_char='c', sequence_row=1, begin_row=2, end_row=3, force=True))

        return outputs


@dxpy.entry_point('run_regenie_step2')
def run_regenie_step2(
        bgen_file: str, bgen_sample: str, bgen_index: str,
        chromosome: str, tarball_prefixes: List[str], samples_include: str,
        covariate_file: Path, pheno_file: Path, pheno_column: str, fit_out_pred: Path, fit_out_loco: Path,
        annotation_file: str, mask_file: str, setlist_file: str, is_binary: bool,
        found_quantitative_covariates: List[str], found_categorical_covariates: List[str],
        ignore_base_covariates: bool) -> Dict[str, Any]:
    """
    A function to run REGENIE step 2 on a single chromosome

    :param bgen_file: The bgen file to run
    :param bgen_sample: The bgen sample file to run
    :param bgen_index: The bgen index file to run
    :param chromosome: The chromosome/chunk to run
    :param tarball_prefixes: The tarball prefixes to run
    :param samples_include: The samples to include
    :param covariate_file: The covariate file to use
    :param pheno_file: The phenotype file to use
    :param pheno_column: The phenotype column to use
    :param fit_out_pred: The prediction file from step 1
    :param fit_out_loco: The LOCO file from step 1
    :param annotation_file: Annotation files
    :param mask_file: Mask files
    :param setlist_file: Setlist files
    :param is_binary: Whether the phenotype is binary
    :param found_quantitative_covariates: Quantitative covariates to use
    :param found_categorical_covariates: Categorical covariates to use
    :param ignore_base_covariates: Whether to ignore base covariates

    :return: The tarball prefix, finished chromosome, and log file
    """

    # we are now working per chunk
    # Get all the files we need
    bgen_file = InputFileHandler(bgen_file).get_file_handle()
    sample_file = InputFileHandler(bgen_sample).get_file_handle()
    samples_include = InputFileHandler(samples_include).get_file_handle()
    bgen_index = InputFileHandler(bgen_index).get_file_handle()
    covariate_file = InputFileHandler(covariate_file).get_file_handle()
    pheno_file = InputFileHandler(pheno_file).get_file_handle()
    fit_out_pred = InputFileHandler(fit_out_pred).get_file_handle()
    InputFileHandler(fit_out_loco).get_file_handle()
    for annotation_file in annotation_file:
        InputFileHandler(annotation_file, download_now=True)
    for mask_file in mask_file:
        InputFileHandler(mask_file, download_now=True)
    for setlist_file in setlist_file:
        InputFileHandler(setlist_file, download_now=True)

    print ("Files in current directory:")
    for path in Path('.').glob('*'):
        LOGGER.info(f'Found file: {path}')

    sample_df = pd.read_csv(sample_file, sep=" ", low_memory=False)
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
        sample_df.to_csv(Path(sample_file), sep=" ", index=False)

    # 4. Run step 2 of regenie
    LOGGER.info("Running REGENIE step 2")
    thread_utility = ThreadUtility(thread_factor=1)

    for tarball_prefix in tarball_prefixes:

        LOGGER.info(f"Running for the mask {tarball_prefix}")
        # print all the files in the current directory
        for path in Path('.').glob('*'):
            LOGGER.info(f'Found file: {path}')
        # if Path(f'{tarball_prefix}.{chromosome}.REGENIE.annotationFile.tsv').exists():
        if Path(f'{tarball_prefix}.{chromosome}.REGENIE.annotationFile.txt').exists():
            thread_utility.launch_job(function=regenie_step_two,
                                      inputs={
                                          "tarball_prefix": tarball_prefix,
                                          "chromosome": chromosome,
                                          "bgen_file": bgen_file,
                                          "bgen_sample": sample_file,
                                          "samples_include": samples_include,
                                          "covariate_file": covariate_file,
                                          "pheno_file": pheno_file,
                                          "pheno_column": pheno_column,
                                          "fit_out_pred": fit_out_pred,
                                          "is_binary": is_binary,
                                          "annotation_file": f'{tarball_prefix}.{chromosome}.REGENIE.annotationFile.txt',
                                          "mask_file": f'{tarball_prefix}.{chromosome}.REGENIE.maskfile.txt',
                                          "setlist_file": f'{tarball_prefix}.{chromosome}.REGENIE.setListFile.txt',
                                          "found_quantitative_covariates": found_quantitative_covariates,
                                          "found_categorical_covariates": found_categorical_covariates,
                                          "ignore_base_covariates": ignore_base_covariates
                                      },
                                      outputs=[
                                          "tarball_prefix",
                                          "chromosome",
                                          "regenie_log",
                                          "regenie_output"
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
            "current_log": exporter.export_files(result["regenie_log"]),
            "regenie_output": exporter.export_files(result["regenie_output"])
        })

    return {
        "output": output
    }


def regenie_step_two(tarball_prefix, chromosome, bgen_file, bgen_sample, samples_include,
                     covariate_file, pheno_file, pheno_column, fit_out_pred, annotation_file,
                     mask_file, setlist_file, found_quantitative_covariates,
                     found_categorical_covariates, is_binary, ignore_base_covariates) -> Tuple[Any, Any, Path, Path]:
    """
    Execution of REGENIE step 2 for a single chromosome

    :param tarball_prefix: The prefix for the tarball files
    :param chromosome: The chromosome/chunk to run
    :param bgen_file: The bgen file to run
    :param bgen_sample: The bgen sample file to run
    :param samples_include: The samples to include
    :param covariate_file: The covariate file to use
    :param pheno_file: The phenotype file to use
    :param pheno_column: The phenotype column to use
    :param fit_out_pred: The prediction file from step 1
    :param annotation_file: The annotation file
    :param mask_file: The mask file
    :param setlist_file: The setlist file
    :param found_quantitative_covariates: The quantitative covariates to use
    :param found_categorical_covariates: The categorical covariates to use
    :param is_binary: Whether the phenotype is binary
    :param ignore_base_covariates: Whether to ignore base covariates

    :return: The tarball prefix, chromosome, and log file
    """

    # set the command executor
    cmd_exec = build_default_command_executor()

    # Note – there is some issue with skato (in --vc-tests flag), so I have changed to skato-acat which works...?
    cmd = f'regenie ' \
          f'--step 2 ' \
          f'--bgen {bgen_file} ' \
          f'--sample {bgen_sample} ' \
          f'--keep {samples_include} ' \
          f'--covarFile {covariate_file} ' \
          f'--phenoFile {pheno_file} ' \
          f'--phenoCol {pheno_column} ' \
          f'--pred {fit_out_pred} ' \
          f'--anno-file {annotation_file} ' \
          f'--mask-def {mask_file} ' \
          f'--set-list {setlist_file} ' \
          f'--aaf-bins 1 ' \
          f'--vc-tests skato-acat,acato-full ' \
          f'--bsize 400 ' \
          f'--threads 1 ' \
          f'--minMAC 1 ' \
          f'--ref-first ' \
          f'--maxCatLevels 110 ' \
          f'--out {tarball_prefix}.{chromosome} '

    cmd += define_covariate_string(found_quantitative_covariates,
                                   found_categorical_covariates,
                                   is_binary,
                                   add_array=False,
                                   ignore_base=ignore_base_covariates)

    regenie_log = Path(f'{tarball_prefix}.{chromosome}.REGENIE_genes.log')
    cmd_exec.run_cmd_on_docker(cmd, stdout_file=regenie_log)

    regenie_output = Path(f'{tarball_prefix}.{chromosome}_{pheno_column}.regenie')

    return tarball_prefix, chromosome, regenie_log, regenie_output
