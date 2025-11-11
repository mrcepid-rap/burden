"""Execute gene-level GLM runs and downstream plotting for burden analyses."""

from pathlib import Path

import pandas as pd
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.linear_model.linear_model import linear_model_null, \
    load_linear_model_genetic_data, run_linear_model
from general_utilities.linear_model.proccess_model_output import process_model_outputs
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner


class GLMRunner(ToolRunner):
    """Run the GLM workflow over all tarball prefixes and available genes."""

    def run_tool(self) -> None:
        """Load input data, fit models, persist results, and generate plots."""

        # 1. Do setup for the linear models.
        # This will load all variants, genes, and phenotypes into memory to allow for parallelization
        # This function returns a class of type LinearModelPack containing info for running GLMs
        self._logger.info("Loading data and running null Linear Model")
        null_model = linear_model_null(self._association_pack.final_covariates,
                                       self._association_pack.pheno_names[0],
                                       self._association_pack.is_binary,
                                       self._association_pack.ignore_base_covariates,
                                       self._association_pack.found_quantitative_covariates,
                                       self._association_pack.found_categorical_covariates)

        # 2. Load the tarballs INTO separate genotypes dictionaries
        self._logger.info("Loading Linear Model genotypes")
        loader_thread_utility = ThreadUtility(threads=self._association_pack.threads,
                                              thread_factor=2,
                                              incrementor=10)

        for tarball_prefix in self._association_pack.tarball_prefixes:
            self._logger.debug("Scheduling genotype load for %s", tarball_prefix)
            loader_thread_utility.launch_job(load_linear_model_genetic_data,
                                             inputs={'tarball_prefix': tarball_prefix,
                                                     'tarball_type': self._association_pack.tarball_type},
                                             outputs=['tarball_prefix', 'genotype_dict'])

        loader_thread_utility.submit_and_monitor()

        # 3. Iterate through every model / gene (in linear_model_pack['genes']) pair and run a GLM
        # We do this directly from the results of the threads above to save storing genotypes in memory
        self._logger.info("Submitting Linear Models to threads")
        runner_thread_utility = ThreadUtility(threads=self._association_pack.threads,
                                              thread_factor=1,
                                              incrementor=500)

        for result in loader_thread_utility:
            for gene in result['genotype_dict'].index.levels[0]:
                runner_thread_utility.launch_job(run_linear_model,
                                                 inputs={'linear_model_pack': null_model,
                                                         'genotype_table': result['genotype_dict'],
                                                         'gene': gene,
                                                         'mask_name': result['tarball_prefix'],
                                                         'is_binary': self._association_pack.is_binary},
                                                 outputs=['gene_dict'])

        runner_thread_utility.submit_and_monitor()

        completed_models = [result['gene_dict'] for result in runner_thread_utility]

        # 4. Annotate unformatted results and print final tabular outputs
        self._logger.info("Processing Linear Model results")
        output_path = Path(f'{self._output_prefix}.genes.glm.stats.tsv')
        self._outputs.extend(process_model_outputs(completed_models,
                                                   output_path,
                                                   self._association_pack.tarball_type,
                                                   self._transcripts_table))

        # 5. Generate Manhattan Plots
        plot_dir = Path(f'{self._output_prefix}_plots')
        plot_dir.mkdir(parents=True, exist_ok=True)
        self._outputs.append(plot_dir)
        glm_table = pd.read_csv(Path(f'{self._output_prefix}.genes.GLM.stats.tsv.gz'), sep='\t')

        mask_values = glm_table['MASK'].value_counts().index
        maf_values = glm_table['MAF'].value_counts().index

        for mask in mask_values:
            for maf in maf_values:
                # To note on the below: I use SYMBOL for the id_column parameter below because ENST is the
                # index and I don't currently have a way to pass the index through to the Plotter methods...
                manhattan_plotter = ManhattanPlotter(self._association_pack.cmd_executor,
                                                     glm_table.query(f'MASK == "{mask}" & MAF == "{maf}"'),
                                                     chrom_column='chrom', pos_column='start',
                                                     alt_column=None,
                                                     id_column='SYMBOL', p_column='p_val_init',
                                                     csq_column='MASK',
                                                     maf_column='cMAC', gene_symbol_column='SYMBOL',
                                                     clumping_distance=1,
                                                     maf_cutoff=30,
                                                     sig_threshold=1E-6)

                manhattan_plotter.plot()[0].rename(plot_dir / f'{mask}.{maf}.genes.GLM.png')
