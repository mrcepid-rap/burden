import csv
import gzip
from pathlib import Path
from typing import List, Tuple

import dxpy
import pandas as pd
from general_utilities.association_resources import define_field_names_from_pandas, bgzip_and_tabix
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import download_bgen_file, LOGGER
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner


class BOLTRunner(ToolRunner):

    def run_tool(self) -> None:

        # Need to pare down the bgen file to samples being tested
        # Have to do this for all chromosomes and all included tarballs, so going to parallelise:

        # 1. First we need to download / prep the BGEN files we want to run through BOLT
        self._logger.info("Processing BGEN files for BOLT run...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=4)

        # The 'poss_chromosomes.txt' has a slightly different format depending on the data-type being used, but
        # generally has a format of <genetics file>\t<fam file>
        possible_chromosomes = Path('poss_chromosomes.txt')
        with open(possible_chromosomes, 'w') as poss_chromosomes:
            for chromosome_chunk in self._association_pack.bgen_dict:
                for tarball_prefix in self._association_pack.tarball_prefixes:
                    if Path(f'{tarball_prefix}.{chromosome_chunk}.BOLT.bgen').exists():
                        thread_utility.launch_job(
                            function=self._process_bolt_bgen_file,
                            inputs={
                                'tarball_prefix': tarball_prefix,
                                'chromosome': chromosome_chunk
                            },
                            outputs=["bgen_output", "sample_file"]
                        )
            thread_utility.submit_and_monitor()

            for result in thread_utility:
                bgen = str(result['bgen_output'])
                sample = str(result['sample_file'])
                poss_chromosomes.write(f'{bgen}\t{sample}\n')

            poss_chromosomes.close()

        # 2. Actually run BOLT
        self._logger.info("Running BOLT...")
        bolt_multithread(poss_chromosomes=possible_chromosomes,
                         bed_file=Path(f"{self._association_pack.genetic_filename}.bed"),
                         bim_file=Path(f"{self._association_pack.genetic_filename}.bim"),
                         fam_file=Path(f"{self._association_pack.genetic_filename}.fam"),
                         low_mac_list=self._association_pack.low_mac_list,
                         final_covariates=self._association_pack.final_covariates,
                         pheno_names=self._association_pack.pheno_names[0],
                         threads=self._association_pack.threads,
                         output_prefix=self._output_prefix,
                         quantitative_covariates=self._association_pack.found_quantitative_covariates,
                         categorical_covariates=self._association_pack.found_categorical_covariates,
                         is_bolt_non_infinite=self._association_pack.is_bolt_non_infinite,
                         ignore_base_covariates=self._association_pack.ignore_base_covariates)

        # 3. Process the outputs
        self._logger.info("Processing BOLT outputs...")
        self._outputs.extend(self._process_bolt_outputs())

    # This handles processing of mask and whole-exome bgen files for input into BOLT
    def _process_bolt_bgen_file(self, tarball_prefix: str, chromosome: str) -> Tuple[Path, Path]:

        # plink2 has implemented a 'fix' for chrX that does not allow samples without sex. This code adds sex back to
        # the bgen sample file so that plink2 will process the data properly. This is the only way I know how to rename
        # variant IDs within a bgen file. For future devs:
        # I code every individual as female. This is because BOLT cannot handle plink2's default coding of ploidy
        # for males as 1. I recognise this is a hack, but for the purposes of this pipeline it is acceptable.
        with Path(f'{tarball_prefix}.{chromosome}.BOLT.sample').open('r') as bgen_sample, \
                Path(f'{tarball_prefix}.{chromosome}.BOLT.fix.sample').open('w') as bgen_fix_sample:

            sample_reader = csv.DictReader(bgen_sample, delimiter=' ')
            sample_writer = csv.DictWriter(bgen_fix_sample, delimiter=' ', fieldnames=sample_reader.fieldnames)
            sample_writer.writeheader()

            # Write sex to 'fix' sample file
            for sample in sample_reader:
                if sample['ID_1'] == '0':
                    sample_writer.writerow(sample)
                else:
                    sample['sex'] = '2'
                    sample_writer.writerow(sample)

        # Do the mask first...
        # We need to modify the bgen file to have an alternate name for IDing masks
        cmd = f'plink2 --threads 4 --bgen {tarball_prefix}.{chromosome}.BOLT.bgen \'ref-first\' ' \
              f'--out {tarball_prefix}.{chromosome} ' \
              f'--make-just-pvar ' \
              f'--sample {tarball_prefix}.{chromosome}.BOLT.fix.sample ' \
              f'--split-par hg38'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        with open(f'{tarball_prefix}.{chromosome}.fixer', 'w') as fix_writer:
            pvar_reader = csv.DictReader(open(f'{tarball_prefix}.{chromosome}.pvar', 'r'), delimiter='\t')
            for variant_id in pvar_reader:
                fix_writer.write(f'{variant_id["ID"]} {variant_id["ID"]}-{tarball_prefix}\n')
            fix_writer.close()

        cmd = f'plink2 --threads 4 --bgen {tarball_prefix}.{chromosome}.BOLT.bgen \'ref-first\' ' \
              f'--sample {tarball_prefix}.{chromosome}.BOLT.fix.sample ' \
              f'--update-name {tarball_prefix}.{chromosome}.fixer ' \
              f'--export bgen-1.2 \'bits=\'8 ' \
              f'--split-par hg38 ' \
              f'--out {tarball_prefix}.{chromosome} '
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        # Make sure the original sample file is being used, otherwise BOLT complains
        sample_file = Path(f'{tarball_prefix}.{chromosome}.BOLT.sample').replace(
            Path(f'{tarball_prefix}.{chromosome}.sample'))

        # set the output file
        bgen_output = Path(f'{tarball_prefix}.{chromosome}.bgen')

        return bgen_output, sample_file

    # This parses the BOLT output file into a usable format for plotting/R
    def _process_bolt_outputs(self) -> List[Path]:

        # First read in the BOLT stats file:
        bolt_table = pd.read_csv(f'{self._output_prefix}.bgen.stats.gz', sep="\t")
        plot_dir = Path(f'{self._output_prefix}_plots/')  # Path to store plots
        plot_dir.mkdir(exist_ok=True)

        # Split the main table into marker and gene tables and remove the larger table
        bolt_table_gene = bolt_table[bolt_table['SNP'].str.contains('ENST')]
        bolt_table_marker = bolt_table[bolt_table['SNP'].str.contains(':')]
        del bolt_table

        # Test what columns we have in the 'SNP' field so we can name them...
        field_names = define_field_names_from_pandas(id_field=bolt_table_gene.iloc[0], default_fields=['ENST'])
        bolt_table_gene[field_names] = bolt_table_gene['SNP'].str.split("-", expand=True)
        bolt_table_gene = bolt_table_gene.drop(columns=['SNP', 'CHR', 'BP', 'ALLELE1', 'ALLELE0', 'GENPOS'])

        # We need to add in an 'AC' column. Pull samples total from the BOLT log file:
        n_bolt = 0
        bolt_log_path = Path(f'{self._output_prefix}.BOLT.log')
        with bolt_log_path.open('r') as bolt_log_file:
            for line in bolt_log_file:
                if 'Total indivs stored in memory:' in line:
                    n_bolt = int(line.strip('Total indivs stored in memory: N = '))
                    break
            bolt_log_file.close()
        # And use them to calculate a MAC
        bolt_table_gene['AC'] = bolt_table_gene['A1FREQ'] * (n_bolt * 2)
        bolt_table_gene['AC'] = bolt_table_gene['AC'].round()

        # Now merge the transcripts table into the gene table to add annotation and the write
        bolt_table_gene = pd.merge(self._transcripts_table, bolt_table_gene, on='ENST', how="left")

        stats_path = Path(f'{self._output_prefix}.genes.BOLT.stats.tsv')
        with stats_path.open('w') as gene_out:
            # Sort by chrom/pos just to be sure...
            bolt_table_gene = bolt_table_gene.sort_values(by=['chrom', 'start', 'end'])
            bolt_table_gene.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')

            # Make Manhattan plots
            for mask in bolt_table_gene['MASK'].value_counts().index:

                for maf in bolt_table_gene['MAF'].value_counts().index:
                    # To note on the below: I use SYMBOL for the id_column parameter below because ENST is the
                    # index and I don't currently have a way to pass the index through to the Plotter methods...
                    manhattan_plotter = ManhattanPlotter(self._association_pack.cmd_executor,
                                                         bolt_table_gene.query(
                                                             f'MASK == "{mask}" & MAF == "{maf}"'),
                                                         chrom_column='chrom', pos_column='start',
                                                         alt_column=None,
                                                         id_column='ENST',
                                                         p_column='P_BOLT_LMM' if self._association_pack.is_bolt_non_infinite else 'P_BOLT_LMM_INF',
                                                         csq_column='MASK',
                                                         maf_column='A1FREQ', gene_symbol_column='SYMBOL',
                                                         clumping_distance=1,
                                                         maf_cutoff=30 / (n_bolt * 2),
                                                         sig_threshold=1E-6)

                    manhattan_plotter.plot()[0].rename(plot_dir / f'{mask}.{maf}.genes.BOLT.png')

        # Define an output(s) array to return
        outputs = [Path(f'{self._output_prefix}.stats.gz'),
                   bolt_log_path,
                   plot_dir]

        # And bgzip and tabix...
        outputs.extend(bgzip_and_tabix(stats_path, skip_row=1, sequence_row=2, begin_row=3, end_row=4))

        # And now process the SNP file (if necessary):
        # Read in the variant index (per-chromosome and mash together)
        if self._association_pack.run_marker_tests:
            variant_index = []
            # Open all chromosome indicies and load them into a list and append them together
            for chromosome_chunk in self._association_pack.bgen_dict:
                variant_index.append(pd.read_csv(f'{chromosome_chunk}.filtered.vep.tsv.gz',
                                                 sep="\t",
                                                 dtype={'SIFT': str, 'POLYPHEN': str}))

            variant_index = pd.concat(variant_index)
            variant_index = variant_index.set_index('varID')

            # For markers, we can use the SNP ID column to get what we need
            bolt_table_marker = bolt_table_marker.rename(columns={'SNP': 'varID', 'A1FREQ': 'BOLT_MAF'})
            bolt_table_marker = bolt_table_marker.drop(columns=['CHR', 'BP', 'ALLELE1', 'ALLELE0', 'GENPOS'])
            bolt_table_marker['BOLT_AC'] = bolt_table_marker['BOLT_MAF'] * (n_bolt * 2)
            bolt_table_marker['BOLT_AC'] = bolt_table_marker['BOLT_AC'].round()
            bolt_table_marker = pd.merge(variant_index, bolt_table_marker, on='varID', how="left")

            marker_tsv = Path(f'{self._output_prefix}.markers.BOLT.stats.tsv')
            with marker_tsv.open('w') as marker_out:

                # Sort by chrom/pos just to be sure...
                bolt_table_marker = bolt_table_marker.sort_values(by=['CHROM', 'POS'])
                bolt_table_marker.to_csv(path_or_buf=marker_out, index=False, sep="\t", na_rep='NA')

            # And bgzip and tabix...
            outputs.extend(bgzip_and_tabix(marker_tsv, skip_row=1, sequence_row=2, begin_row=3))

        return outputs


def bolt_multithread(poss_chromosomes: Path, bed_file: Path, bim_file: Path, fam_file: Path, low_mac_list: Path,
                     final_covariates: Path,
                     pheno_names: str, threads: int, output_prefix: str, quantitative_covariates: List,
                     categorical_covariates: List,
                     is_bolt_non_infinite: bool, ignore_base_covariates: bool) -> None:
    """
    Multithread BOLT runner for WGS data.

    :param poss_chromosomes: Path to the file containing the list of BGEN files and their sample files.
    :param bed_file: Path to the plink .bed file of array data
    :param bim_file: Path to the plink .bim file of array data
    :param fam_file: Path to the plink .fam file of array data
    :param low_mac_list: Path to the list of low MAC markers to exclude
    :param final_covariates: Final covariate file to use
    :param pheno_names: List of phenotype names to test
    :param threads: Number of threads to use
    :param output_prefix: Output prefix to use
    :param quantitative_covariates: List of quantitative covariates found in the covariate file
    :param categorical_covariates: List of categorical covariates found in the covariate file
    :param is_bolt_non_infinite: Whether to run BOLT non-infinitesimal
    :param ignore_base_covariates: Whether to ignore the base covariates
    :return: None
    """

    LOGGER.info("Starting BOLT run using subjobs...")

    # read in all the possible chromosomes/chunks that we can run
    poss_chromosomes = pd.read_csv(poss_chromosomes, sep='\t', header=None)

    # set the subjob launcher class
    subjob_launcher = joblauncher_factory()

    # set the file exporter for each chunk
    exporter = ExportFileHandler()
    bed_file = exporter.export_files(bed_file)
    bim_file = exporter.export_files(bim_file)
    fam_file = exporter.export_files(fam_file)
    low_mac_list = exporter.export_files(low_mac_list)
    final_covariates = exporter.export_files(final_covariates)

    for idx, row in poss_chromosomes.iterrows():
        # write each row as a new file, so that we can launch a subjob for each
        chromosome_file = Path(f'chromosome_{idx}.txt')
        chromosome_file.write_text(f'{row.iloc[0]}\t{row.iloc[1]}\n')
        chromosome_file = exporter.export_files(chromosome_file)
        bgen_file = exporter.export_files(row.iloc[0])
        sample_file = exporter.export_files(row.iloc[1])

        print(bgen_file)
        print(sample_file)

        return

        subjob_launcher.launch_job(
            function=run_bolt,
            inputs={
                'possible_chromosomes': chromosome_file,
                'bgen_file': bgen_file,
                'sample_file': sample_file,
                'bed_file': bed_file,
                'bim_file': bim_file,
                'fam_file': fam_file,
                'low_mac_list': low_mac_list,
                'final_covariates': final_covariates,
                'pheno_names': pheno_names,
                'threads': threads,
                'output_prefix': f'{output_prefix}_part{idx}',
                'found_quantitative_covariates': quantitative_covariates,
                'found_categorical_covariates': categorical_covariates,
                "is_bolt_non_infinite": is_bolt_non_infinite,
                'ignore_base_covariates': ignore_base_covariates
            },
            outputs=['output'],
            name=f"BOLT chunk {idx}"
        )
    subjob_launcher.submit_and_monitor()

    LOGGER.info("BOLT subjobs finished successfully.")

    # process results
    for subjob_output in subjob_launcher:
        output = subjob_output['output']
        InputFileHandler(output["bolt_log"], download_now=True)
        InputFileHandler(output["statsfile"], download_now=True)
        InputFileHandler(output["statsfile_bgen_snps"], download_now=True)

    # get *.stats.gz (but not bgen stats)
    stats_parts = sorted([
        f for f in Path().glob(f"{output_prefix}_part*.stats.gz")
        if not f.name.endswith(".bgen.stats.gz")
    ])
    # get *.bgen.stats.gz
    bgen_stats_parts = sorted(Path().glob(f"{output_prefix}_part*.bgen.stats.gz"))
    # get *.BOLT.log
    log_parts = sorted(Path().glob(f"{output_prefix}_part*.BOLT.log"))

    # Concatenate .stats.gz
    with gzip.open(f"{output_prefix}.stats.gz", "wb") as stats:
        for i, part_file in enumerate(stats_parts):
            with gzip.open(part_file, "rb") as fin:
                if i > 0:
                    fin.readline()  # skip header
                stats.write(fin.read())
    # Concatenate .bgen.stats.gz
    with gzip.open(f"{output_prefix}.bgen.stats.gz", "wb") as bgen_stats:
        for i, part_file in enumerate(bgen_stats_parts):
            with gzip.open(part_file, "rb") as fin:
                if i > 0:
                    fin.readline()
                bgen_stats.write(fin.read())
    # Concatenate .BOLT.log
    with open(f"{output_prefix}.BOLT.log", "w") as log:
        for part_file in log_parts:
            with open(part_file, "r") as fin:
                log.write(fin.read())
                log.write("\n")  # ensure spacing between parts


# Run rare variant association testing using BOLT
@dxpy.entry_point('run_bolt')
def run_bolt(possible_chromosomes: str, bgen_file: str, sample_file: str,
             bed_file: str, bim_file: str, fam_file: str, low_mac_list: str, final_covariates: str,
             pheno_names: str, threads: int, output_prefix: str, found_quantitative_covariates: List[str],
             found_categorical_covariates: List[str], is_bolt_non_infinite: bool, ignore_base_covariates: bool) -> dict:
    """
    Run BOLT on the WES data using the covariates and genetic data provided.

    :param possible_chromosomes: Path to the file containing the list of BGEN files and their sample files.
    :param bgen_file: Path to the BGEN file to use
    :param sample_file: Path to the sample file to use
    :param bed_file: Path to the plink .bed file of array data
    :param bim_file: Path to the plink .bim file of array data
    :param fam_file: Path to the plink .fam file of array data
    :param low_mac_list:
    :param final_covariates: Final covariate file to use
    :param pheno_names: List of phenotype names to test
    :param threads: Number of threads to use
    :param output_prefix: Output prefix to use
    :param found_quantitative_covariates: List of quantitative covariates found in the covariate file
    :param found_categorical_covariates: List of categorical covariates found in the covariate file
    :param is_bolt_non_infinite: Whether to run BOLT in non-infinitesimal mode
    :param ignore_base_covariates: Whether to ignore the base covariates

    :return: None
    """

    # Download all the files we need locally
    possible_chromosomes = InputFileHandler(possible_chromosomes).get_file_handle()
    InputFileHandler(bgen_file).get_file_handle()
    InputFileHandler(sample_file).get_file_handle()
    bed_file = InputFileHandler(bed_file).get_file_handle()
    bim_file = InputFileHandler(bim_file).get_file_handle()
    fam_file = InputFileHandler(fam_file).get_file_handle()
    low_mac_list = InputFileHandler(low_mac_list).get_file_handle()
    final_covariates = InputFileHandler(final_covariates).get_file_handle()

    # set the command executor
    cmd_exec = build_default_command_executor()

    # set the outputs
    bolt_log = Path(f'{output_prefix}.BOLT.log')
    statsfile = Path(f'{output_prefix}.stats.gz')
    statsfile_bgen_snps = Path(f'{output_prefix}.bgen.stats.gz')

    # See the README.md for more information on these parameters
    # REMEMBER: The geneticMapFile is for the bfile, not the WES data!
    # REMEMBER we do --noBgenIDcheck because the genetic data is filtered to the covariate file, the bgens are NOT
    cmd = f'bolt ' + \
          f'--bed={bed_file} ' \
          f'--bim={bim_file} ' \
          f'--fam={fam_file} ' \
          f'--exclude={low_mac_list} ' \
          f'--phenoFile={final_covariates} ' \
          f'--phenoCol={pheno_names} ' \
          f'--covarFile={final_covariates} ' \
          f'--covarMaxLevels=110 ' \
          f'--LDscoresFile=/home/app/BOLT-LMM_v2.4.1/tables/LDSCORE.1000G_EUR.tab.gz ' \
          f'--geneticMapFile=/home/app/BOLT-LMM_v2.4.1/tables/genetic_map_hg19_withX.txt.gz ' \
          f'--numThreads={threads} ' \
          f'--statsFile={statsfile} ' \
          f'--verboseStats ' \
          f'--bgenSampleFileList={possible_chromosomes} ' \
          f'--noBgenIDcheck ' \
          f'--LDscoresMatchBp ' \
          f'--statsFileBgenSnps={statsfile_bgen_snps} '

    if is_bolt_non_infinite:
        cmd += '--lmmForceNonInf '
    else:
        cmd += '--lmmInfOnly '

    if not ignore_base_covariates:
        cmd += f'--covarCol=sex ' \
               f'--covarCol=wes_batch ' \
               f'--qCovarCol=age ' \
               f'--qCovarCol=age_squared ' \
               f'--qCovarCol=PC{{1:10}} '

    if len(found_quantitative_covariates) > 0:
        for covar in found_quantitative_covariates:
            cmd += f'--qCovarCol={covar} '
    if len(found_categorical_covariates) > 0:
        for covar in found_categorical_covariates:
            cmd += f'--covarCol={covar} '

    cmd_exec.run_cmd_on_docker(cmd, stdout_file=bolt_log)

    # check that all the output files were created
    for file in [bolt_log, statsfile, statsfile_bgen_snps]:
        if not file.exists():
            raise dxpy.AppError(f'BOLT failed to produce expected output file: {file}')

    exporter = ExportFileHandler()
    output = {
        'bolt_log': exporter.export_files(bolt_log),
        'statsfile': exporter.export_files(statsfile),
        'statsfile_bgen_snps': exporter.export_files(statsfile_bgen_snps),
    }

    return output
