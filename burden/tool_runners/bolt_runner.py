import csv
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from general_utilities.association_resources import define_field_names_from_pandas, bgzip_and_tabix
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
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
        # print current directory contents for debugging
        for file in Path('.').iterdir():
            self._logger.debug(f'Found file: {file.name}')
        with open(possible_chromosomes, 'w') as poss_chromosomes:
            for chromosome_chunk in self._association_pack.bgen_dict:
                print(chromosome_chunk)
                for tarball_prefix in self._association_pack.tarball_prefixes:
                    if not Path(f'{tarball_prefix}.{chromosome_chunk}.BOLT.bgen').exists():
                        print(f'{chromosome_chunk} does not exist')
                        continue
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
                bgen, sample = result.values()
                poss_chromosomes.write(f'{bgen} '
                                       f'{sample}\n')

            poss_chromosomes.close()

        # 2. Actually run BOLT
        self._logger.info("Running BOLT...")
        self._run_bolt(possible_chromosomes)

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

    # Run rare variant association testing using BOLT
    def _run_bolt(self, possible_chromosomes: Path) -> None:
        """
        Run BOLT on the WES data using the covariates and genetic data provided.

        :param possible_chromosomes: Path to the file containing the list of BGEN files and their sample files.
        :return: None
        """

        # See the README.md for more information on these parameters
        # REMEMBER: The geneticMapFile is for the bfile, not the WES data!
        # REMEMBER we do --noBgenIDcheck because the genetic data is filtered to the covariate file, the bgens are NOT
        cmd = f'bolt ' + \
              f'--bed={self._association_pack.genetic_filename}.bed ' \
              f'--bim={self._association_pack.genetic_filename}.bim ' \
              f'--fam={self._association_pack.genetic_filename}.fam ' \
              f'--exclude={self._association_pack.low_mac_list} ' \
              f'--phenoFile={self._association_pack.final_covariates} ' \
              f'--phenoCol={self._association_pack.pheno_names[0]} ' \
              f'--covarFile={self._association_pack.final_covariates} ' \
              f'--covarMaxLevels=110 ' \
              f'--LDscoresFile=/home/app/BOLT-LMM_v2.4.1/tables/LDSCORE.1000G_EUR.tab.gz ' \
              f'--geneticMapFile=/home/app/BOLT-LMM_v2.4.1/tables/genetic_map_hg19_withX.txt.gz ' \
              f'--numThreads={self._association_pack.threads} ' \
              f'--statsFile={self._output_prefix}.stats.gz ' \
              f'--verboseStats ' \
              f'--bgenSampleFileList={possible_chromosomes.name} ' \
              f'--noBgenIDcheck ' \
              f'--LDscoresMatchBp ' \
              f'--statsFileBgenSnps={self._output_prefix}.bgen.stats.gz '

        if self._association_pack.is_bolt_non_infinite:
            cmd += '--lmmForceNonInf '
        else:
            cmd += '--lmmInfOnly '

        if not self._association_pack.ignore_base_covariates:
            cmd += f'--covarCol=sex ' \
                   f'--covarCol=batch ' \
                   f'--qCovarCol=age ' \
                   f'--qCovarCol=age_squared ' \
                   f'--qCovarCol=PC{{1:10}} '

        if len(self._association_pack.found_quantitative_covariates) > 0:
            for covar in self._association_pack.found_quantitative_covariates:
                cmd += f'--qCovarCol={covar} '
        if len(self._association_pack.found_categorical_covariates) > 0:
            for covar in self._association_pack.found_categorical_covariates:
                cmd += f'--covarCol={covar} '
        bolt_log = Path(f'{self._output_prefix}.BOLT.log')

        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=bolt_log)

    # This parses the BOLT output file into a usable format for plotting/R
    def _process_bolt_outputs(self) -> List[Path]:

        pd.set_option('display.max_rows', None)

        # print all files in the current directory for debugging
        from pathlib import Path
        current_dir = Path('.')
        print("Current directory contents:")
        for file in current_dir.iterdir():
            print(file.name)

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
                                                         bolt_table_gene.query(f'MASK == "{mask}" & MAF == "{maf}"'),
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
            # # Open all chromosome indicies and load them into a list and append them together
            # for chromosome_chunk in self._association_pack.bgen_dict:
            #     variant_index.append(pd.read_csv(f'{chromosome_chunk}.filtered.vep.tsv.gz',
            #                                      sep="\t",
            #                                      dtype={'SIFT': str, 'POLYPHEN': str}))
            #
            # Open all chromosome indices (or single index if one chromosome processing) and load them into a list
            # and append them together
            for chromosome_name, chromosome_info in self._association_pack.bgen_dict.items():
                print(chromosome_name)
                print(chromosome_info)
                local_vep = chromosome_info['vep'].get_file_handle()
                variant_index.append(
                    pd.read_csv(local_vep, sep="\t", dtype={'SIFT': str, 'POLYPHEN': str, }))

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
