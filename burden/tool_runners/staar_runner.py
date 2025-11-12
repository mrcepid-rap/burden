from pathlib import Path
from typing import List

import pandas as pd
from general_utilities.job_management.command_executor import build_default_command_executor, DockerMount
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.linear_model import linear_model
from general_utilities.linear_model.proccess_model_output import process_model_outputs
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner

# import the models
staar_null_script = Path(__file__).parent / 'R_resources' / 'runSTAAR_Null.R'
staar_genes_script = Path(__file__).parent / 'R_resources' / 'runSTAAR_Genes.R'


class STAARRunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Run the STAAR NULL model
        self._logger.info("Running STAAR Null Model...")
        if not Path('phenotype.STAAR_null.rds').exists():
            null_model = linear_model.linear_model_null(
                phenofile=self._association_pack.final_covariates,
                phenotype=self._association_pack.pheno_names[0],
                is_binary=self._association_pack.is_binary,
                ignore_base=self._association_pack.ignore_base_covariates,
                found_categorical_covariates=self._association_pack.found_categorical_covariates,
                found_quantitative_covariates=self._association_pack.found_quantitative_covariates
            )
        else:
            self._logger.info("STAAR null model already exists, skipping...")

        # 2. Run the actual per-gene association tests
        self._logger.info("Running STAAR masks * chromosomes...")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=1)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                for chromosome in self._association_pack.bgen_dict:
                    continue


                    #
                    #     thread_utility.launch_job(self.staar_genes,
                    #                               tarball_prefix=tarball_prefix,
                    #                               chromosome=chromosome,
                    #                               phenoname=phenoname,
                    #                               has_gene_info=False,
                    #                               r_script=staar_genes_script)
                    # else:
                    #     # If the gene matrix file does not exist, we can skip this chromosome
                    #     self._logger.warning(f'No gene matrix file found for {tarball_prefix}.{chromosome}. '
                    #                          f'STAAR analysis will not work.')
                    #     continue

        # 3. Print a preliminary STAAR output
        self._logger.info("Finalising STAAR outputs...")
        completed_staar_files = []
        # And gather the resulting futures
        for result in thread_utility:
            tarball_prefix, finished_chromosome, phenoname = result
            completed_staar_files.append(f'{tarball_prefix}.{phenoname}.{finished_chromosome}.STAAR_results.tsv')

        # 4. Annotate and print final STAAR output
        # self._outputs.extend(process_staar_outputs(completed_staar_files, self._output_prefix))

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

# def multithread_staar_burden():