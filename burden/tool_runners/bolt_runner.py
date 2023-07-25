import csv
import pandas as pd

from pathlib import Path
from typing import List

from burden.tool_runners.tool_runner import ToolRunner
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.association_resources import get_chromosomes, process_bgen_file, \
    build_transcript_table, define_field_names_from_pandas, bgzip_and_tabix


class BOLTRunner(ToolRunner):

    def run_tool(self) -> None:

        # Need to pare down the bgen file to samples being tested
        # Have to do this for all chromosomes and all included tarballs, so going to parallelise:

        # 1. First we need to download / prep the BGEN files we want to run through BOLT
        self._logger.info("Processing BGEN files for BOLT run...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A BOLT thread failed',
                                       incrementor=10,
                                       thread_factor=4)

        # The 'poss_chromosomes.txt' has a slightly different format depending on the data-type being used, but
        # generally has a format of <genetics file>\t<fam file>
        with open('poss_chromosomes.txt', 'w') as poss_chromosomes:
            for chromosome in get_chromosomes():
                for tarball_prefix in self._association_pack.tarball_prefixes:
                    if Path(f'{tarball_prefix}.{chromosome}.BOLT.bgen').exists():
                        poss_chromosomes.write(f'/test/{tarball_prefix}.{chromosome}.bgen '
                                               f'/test/{tarball_prefix}.{chromosome}.sample\n')
                        thread_utility.launch_job(class_type=self._process_bolt_bgen_file,
                                                  tarball_prefix=tarball_prefix,
                                                  chromosome=chromosome)

                if self._association_pack.run_marker_tests:
                    poss_chromosomes.write(f'/test/{chromosome}.markers.bgen '
                                           f'/test/{chromosome}.markers.sample\n')
                    # This makes use of a utility class from AssociationResources since bgen filtering/processing is
                    # IDENTICAL to that done for SAIGE. Do not want to duplicate code!
                    thread_utility.launch_job(class_type=process_bgen_file,
                                              chrom_bgen_index=self._association_pack.bgen_dict[chromosome],
                                              chromosome=chromosome)

            poss_chromosomes.close()
            thread_utility.collect_futures()

        # 2. Actually run BOLT
        self._logger.info("Running BOLT...")
        self._run_bolt()

        # 3. Process the outputs
        self._logger.info("Processing BOLT outputs...")
        self._outputs.extend(self._process_bolt_outputs())

    # This handles processing of mask and whole-exome bgen files for input into BOLT
    def _process_bolt_bgen_file(self, tarball_prefix: str, chromosome: str) -> None:

        # Do the mask first...
        # We need to modify the bgen file to have an alternate name for IDing masks
        cmd = f'plink2 --threads 4 --bgen /test/{tarball_prefix}.{chromosome}.BOLT.bgen \'ref-last\' ' \
                    f'--out /test/{tarball_prefix}.{chromosome} ' \
                    f'--make-just-pvar'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

        with open(f'{tarball_prefix}.{chromosome}.fixer', 'w') as fix_writer:
            pvar_reader = csv.DictReader(open(f'{tarball_prefix}.{chromosome}.pvar', 'r'), delimiter='\t')
            for variant_id in pvar_reader:
                fix_writer.write(f'{variant_id["ID"]} {variant_id["ID"]}-{tarball_prefix}\n')
            fix_writer.close()

        cmd = f'plink2 --threads 4 --bgen /test/{tarball_prefix}.{chromosome}.BOLT.bgen \'ref-last\' ' \
              f'--sample /test/{tarball_prefix}.{chromosome}.BOLT.sample ' \
              f'--update-name /test/{tarball_prefix}.{chromosome}.fixer ' \
              f'--export bgen-1.2 \'bits=\'8 ' \
              f'--out /test/{tarball_prefix}.{chromosome} ' \
              f'--keep-fam /test/SAMPLES_Include.txt'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)

    # Run rare variant association testing using BOLT
    def _run_bolt(self) -> None:

        # See the README.md for more information on these parameters
        # REMEMBER: The geneticMapFile is for the bfile, not the WES data!
        cmd = f'bolt ' + \
                f'--bfile=/test/genetics/UKBB_470K_Autosomes_QCd_WBA ' \
                f'--exclude=/test/genetics/UKBB_470K_Autosomes_QCd.low_MAC.snplist ' \
                f'--phenoFile=/test/phenotypes_covariates.formatted.txt ' \
                f'--phenoCol={self._association_pack.pheno_names[0]} ' \
                f'--covarFile=/test/phenotypes_covariates.formatted.txt ' \
                f'--covarCol=sex ' \
                f'--covarCol=wes_batch ' \
                f'--qCovarCol=age ' \
                f'--qCovarCol=age_squared ' \
                f'--qCovarCol=PC{{1:10}} ' \
                f'--covarMaxLevels=110 ' \
                f'--LDscoresFile=BOLT-LMM_v2.4.1/tables/LDSCORE.1000G_EUR.tab.gz ' \
                f'--geneticMapFile=BOLT-LMM_v2.4.1/tables/genetic_map_hg19_withX.txt.gz ' \
                f'--numThreads={self._association_pack.threads} ' \
                f'--statsFile=/test/{self._output_prefix}.stats.gz ' \
                f'--verboseStats '

        cmd += f'--bgenSampleFileList=/test/poss_chromosomes.txt ' \
               f'--statsFileBgenSnps=/test/{self._output_prefix}.bgen.stats.gz'

        if self._association_pack.is_bolt_non_infinite:
            cmd += ' --lmmForceNonInf'
        else:
            cmd += ' --lmmInfOnly'

        if len(self._association_pack.found_quantitative_covariates) > 0:
            for covar in self._association_pack.found_quantitative_covariates:
                cmd += f' --qCovarCol={covar} '
        if len(self._association_pack.found_categorical_covariates) > 0:
            for covar in self._association_pack.found_categorical_covariates:
                cmd += f' --covarCol={covar} '
        bolt_log = Path(f'{self._output_prefix}.BOLT.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=bolt_log)

    # This parses the BOLT output file into a useable format for plotting/R
    def _process_bolt_outputs(self) -> List[Path]:

        # First read in the BOLT stats file:
        bolt_table = pd.read_csv(f'{self._output_prefix}.bgen.stats.gz', sep="\t")

        # Split the main table into marker and gene tables and remove the larger table
        bolt_table_gene = bolt_table[bolt_table['SNP'].str.contains('ENST')]
        bolt_table_marker = bolt_table[bolt_table['SNP'].str.contains(':')]
        del bolt_table

        # Now process the gene table into a useable format:
        # First read in the transcripts file
        transcripts_table = build_transcript_table()

        # Test what columns we have in the 'SNP' field so we can name them...
        field_names = define_field_names_from_pandas(bolt_table_gene.iloc[0])
        bolt_table_gene[field_names] = bolt_table_gene['SNP'].str.split("-", expand=True)
        bolt_table_gene = bolt_table_gene.drop(columns=['SNP', 'CHR', 'BP', 'ALLELE1', 'ALLELE0', 'GENPOS'])

        # We need to add in an 'AC' column. Pull samples total from the BOLT log file:
        n_bolt = 0
        bolt_log_path = Path(f'{self._output_prefix}.BOLT.log')
        with bolt_log_path.open('r') as bolt_log_file:
            for line in bolt_log_file:
                if 'samples (Nbgen):' in line:
                    n_bolt = int(line.strip('samples (Nbgen): '))
                    break
            bolt_log_file.close()
        # And use them to calculate a MAC
        bolt_table_gene['AC'] = bolt_table_gene['A1FREQ'] * (n_bolt*2)
        bolt_table_gene['AC'] = bolt_table_gene['AC'].round()

        # Now merge the transcripts table into the gene table to add annotation and the write
        bolt_table_gene = pd.merge(transcripts_table, bolt_table_gene, on='ENST', how="left")

        stats_path = Path(f'{self._output_prefix}.genes.BOLT.stats.tsv')
        with stats_path.open('w') as gene_out:
            # Sort by chrom/pos just to be sure...
            bolt_table_gene = bolt_table_gene.sort_values(by=['chrom', 'start', 'end'])
            bolt_table_gene.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')

        # Define an output(s) array to return
        outputs = [Path(f'{self._output_prefix}.stats.gz'),
                   bolt_log_path]

        # And bgzip and tabix...
        outputs.extend(bgzip_and_tabix(stats_path, skip_row=1, sequence_row=2, begin_row=3, end_row=4))

        # And now process the SNP file (if necessary):
        # Read in the variant index (per-chromosome and mash together)
        if self._association_pack.run_marker_tests:
            variant_index = []
            # Open all chromosome indicies and load them into a list and append them together
            for chromosome in get_chromosomes():
                variant_index.append(pd.read_csv(f'filtered_bgen/{chromosome}.filtered.vep.tsv.gz',
                                                 sep="\t",
                                                 dtype={'SIFT': str, 'POLYPHEN': str}))

            variant_index = pd.concat(variant_index)
            variant_index = variant_index.set_index('varID')

            # For markers, we can use the SNP ID column to get what we need
            bolt_table_marker = bolt_table_marker.rename(columns={'SNP': 'varID', 'A1FREQ': 'BOLT_MAF'})
            bolt_table_marker = bolt_table_marker.drop(columns=['CHR', 'BP', 'ALLELE1', 'ALLELE0', 'GENPOS'])
            bolt_table_marker['BOLT_AC'] = bolt_table_marker['BOLT_MAF'] * (n_bolt*2)
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
