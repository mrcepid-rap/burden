import csv
import pandas as pd
from pathlib import Path

from burden.tool_runners.tool_runner import ToolRunner
from general_utilities.association_resources import bgzip_and_tabix
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.linear_model import linear_model
from general_utilities.linear_model.linear_model import LinearModelResult
from general_utilities.linear_model.proccess_model_output import process_model_outputs
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter


class GLMRunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Do setup for the linear models.
        # This will load all variants, genes, and phenotypes into memory to allow for parallelization
        # This function returns a class of type LinearModelPack containing info for running GLMs
        self._logger.info("Loading data and running null Linear Model")
        null_model = linear_model.linear_model_null(self._association_pack.pheno_names[0],
                                                    self._association_pack.is_binary,
                                                    self._association_pack.ignore_base_covariates,
                                                    self._association_pack.found_quantitative_covariates,
                                                    self._association_pack.found_categorical_covariates)

        # 2. Load the tarballs INTO separate genotypes dictionaries
        self._logger.info("Loading Linear Model genotypes")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=10,
                                       thread_factor=2)

        for tarball_prefix in self._association_pack.tarball_prefixes:
            thread_utility.launch_job(linear_model.load_tarball_linear_model,
                                      tarball_prefix=tarball_prefix,
                                      is_snp_tar=False,
                                      is_gene_tar=False)

        genotype_packs = {}
        for result in thread_utility:
            tarball_prefix, genotype_dict = result
            genotype_packs[tarball_prefix] = genotype_dict

        # 3. Iterate through every model / gene (in linear_model_pack['genes']) pair and run a GLM
        self._logger.info("Submitting Linear Models to threads")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       error_message='A GLM thread failed',
                                       incrementor=500,
                                       thread_factor=1)

        for model in genotype_packs:
            for gene in genotype_packs[model].index.levels[0]:  # level[0] in this DataFrame is ENST
                thread_utility.launch_job(linear_model.run_linear_model,
                                          linear_model_pack=null_model,
                                          genotype_table=genotype_packs[model],
                                          gene=gene,
                                          mask_name=model,
                                          is_binary=self._association_pack.is_binary)

        # As futures finish, write unformatted results:
        fieldnames = ['ENST', 'maskname', 'pheno_name', 'p_val_init', 'n_car', 'cMAC', 'n_model',
                      'p_val_full', 'effect', 'std_err']
        # Binary traits get an additional set of fields to describe the confusion matrix.
        if self._association_pack.is_binary:
            fieldnames.extend(['n_noncar_affected', 'n_noncar_unaffected', 'n_car_affected', 'n_car_unaffected'])

        lm_stats_path = Path(f'{self._output_prefix}.lm_stats.tmp')
        with lm_stats_path.open('w') as lm_stats_file:
            lm_stats_csv = csv.DictWriter(lm_stats_file,
                                          delimiter="\t",
                                          fieldnames=fieldnames,
                                          extrasaction='ignore')

            lm_stats_csv.writeheader()
            for result in thread_utility:
                finished_gene: LinearModelResult = result
                lm_stats_csv.writerow(finished_gene.todict())

        # 4. Annotate unformatted results and print final tabular outputs
        self._logger.info("Annotating Linear Model results")
        self._outputs.extend(process_linear_model_outputs(self._output_prefix))

        # 5. Generate Manhattan Plots
        plot_dir = Path(f'{self._output_prefix}_plots/')  # Path to store plots
        plot_dir.mkdir()
        self._outputs.extend(plot_dir)
        glm_table = pd.read_csv(Path(f'{self._output_prefix}.genes.glm.stats.tsv.gz'), sep='\t', )

        for mask in glm_table['MASK'].value_counts().index:

            for maf in glm_table['MAF'].value_counts().index:
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




