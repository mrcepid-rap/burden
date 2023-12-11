from pathlib import Path

import pandas as pd

from burden.tool_runners.tool_runner import ToolRunner
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.association_resources import get_chromosomes
from general_utilities.linear_model.proccess_model_output import process_staar_outputs
from general_utilities.linear_model.staar_model import staar_null, staar_genes
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter


class STAARRunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Run the STAAR NULL model
        self._logger.info("Running STAAR Null Model...")
        staar_null(phenoname=self._association_pack.pheno_names[0],
                   is_binary=self._association_pack.is_binary,
                   sex=self._association_pack.sex,
                   ignore_base=self._association_pack.ignore_base_covariates,
                   found_quantitative_covariates=self._association_pack.found_quantitative_covariates,
                   found_categorical_covariates=self._association_pack.found_categorical_covariates)

        # 2. Run the actual per-gene association tests
        self._logger.info("Running STAAR masks * chromosomes...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A STAAR thread failed',
                                       incrementor=10,
                                       thread_factor=1)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                for chromosome in get_chromosomes():
                    if Path(f'{tarball_prefix}.{chromosome}.STAAR.matrix.rds'):
                        thread_utility.launch_job(staar_genes,
                                                  tarball_prefix=tarball_prefix,
                                                  chromosome=chromosome,
                                                  phenoname=phenoname,
                                                  has_gene_info=False)

        # 3. Print a preliminary STAAR output
        self._logger.info("Finalising STAAR outputs...")
        completed_staar_files = []
        # And gather the resulting futures
        for result in thread_utility:
            tarball_prefix, finished_chromosome, phenoname = result
            completed_staar_files.append(f'{tarball_prefix}.{phenoname}.{finished_chromosome}.STAAR_results.tsv')

        # 4. Annotate and print final STAAR output
        self._outputs.extend(process_staar_outputs(completed_staar_files, self._output_prefix))

        # 5. Make Manhattan plots
        plot_dir = Path(f'{self._output_prefix}_plots/')  # Path to store plots
        plot_dir.mkdir()
        self._outputs.append(plot_dir)

        staar_table_gene = pd.read_csv(f'{self._output_prefix}.genes.STAAR.stats.tsv', sep='\t')

        for mask in staar_table_gene['MASK'].value_counts().index:

            for maf in staar_table_gene['MAF'].value_counts().index:
                # To note on the below: I use SYMBOL for the id_column parameter below because ENST is the
                # index and I don't currently have a way to pass the index through to the Plotter methods...
                manhattan_plotter = ManhattanPlotter(self._association_pack.cmd_executor,
                                                     staar_table_gene.query(f'MASK == "{mask}" & MAF == "{maf}"'),
                                                     chrom_column='chrom', pos_column='start',
                                                     alt_column=None,
                                                     id_column='ENST',
                                                     p_column='staar.O.p',
                                                     csq_column='MASK',
                                                     maf_column='cMAC', gene_symbol_column='SYMBOL',
                                                     clumping_distance=1,
                                                     maf_cutoff=30,
                                                     sig_threshold=1E-6)

                manhattan_plotter.plot()[0].rename(plot_dir / f'{mask}.{maf}.genes.STAAR.png')
