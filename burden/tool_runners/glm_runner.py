"""Execute gene-level GLM runs and downstream plotting for burden analyses."""

import pickle
from pathlib import Path
from typing import List, Any, Dict

import dxpy
import pandas as pd
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.linear_model.linear_model import linear_model_null, \
    load_linear_model_genetic_data, run_linear_model
from general_utilities.linear_model.proccess_model_output import process_model_outputs
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner


class GLMRunner(ToolRunner):
    """Run the GLM workflow over all tarball prefixes and available genes."""

    def run_tool(self) -> None:
        """Load input data, fit models, persist results, and generate plots.

        This method coordinates the entire GLM workflow, including:
        1. Running the null linear model.
        2. Loading genetic data for each tarball prefix.
        3. Launching parallel subjobs to run gene-level GLMs.
        4. Processing and annotating the results.
        5. Generating Manhattan plots.
        """

        # 1. Do setup for the linear models & get the null model
        self._logger.info("Loading data and running null Linear Model")
        null_model = linear_model_null(self._association_pack.final_covariates,
                                       self._association_pack.pheno_names[0],
                                       self._association_pack.is_binary,
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
        genotype_results = list(loader_thread_utility)

        # 3. Iterate through every model / gene pair and run a GLM in parallel subjobs
        self._logger.info("Submitting Linear Models to subjobs")
        completed_models = self._multithread_glm_subjobs(null_model, genotype_results)

        # 4. Annotate unformatted results and print final tabular outputs
        self._logger.info("Processing Linear Model results")
        output_path = Path(f'{self._output_prefix}.genes.GLM.stats.tsv')
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

    def _multithread_glm_subjobs(self, null_model: Any, genotype_results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Launch parallel subjobs to run gene-level GLMs for each tarball prefix.

        This method serializes the null model and genetic data, exports them to DNAnexus,
        and launches subjobs. It then collects and aggregates the results from all subjobs.

        :param null_model: The fitted null linear model object.
        :param genotype_results: A list of dictionaries, where each dictionary contains
                                 'tarball_prefix' and 'genotype_dict' (pandas DataFrame)
                                 for a specific tarball.
        :return: A list of dictionaries, where each dictionary represents the results for a gene.
        """
        # Set up the job launcher
        launcher = joblauncher_factory(download_on_complete=True)
        exporter = ExportFileHandler(delete_on_upload=False)

        # First, serialize the null model
        null_model_path = Path("null_model.pkl")
        with null_model_path.open('wb') as null_writer:
            pickle.dump(null_model, null_writer)
        dnanexus_null_model = exporter.export_files(null_model_path)

        # Launch a subjob for each tarball prefix
        for result in genotype_results:
            # Serialize the genotype dictionary for this tarball prefix
            tarball_prefix = result['tarball_prefix']
            genotype_dict = result['genotype_dict']
            genotype_path = Path(f"{tarball_prefix}.genotypes.pkl")
            genotype_dict.to_pickle(genotype_path)
            dnanexus_genotype_dict = exporter.export_files(genotype_path)

            launcher.launch_job(
                function=run_glm_gene_subjob,
                inputs={
                    'null_model_dxfile': dnanexus_null_model,
                    'genotype_dict_dxfile': dnanexus_genotype_dict,
                    'tarball_prefix': tarball_prefix,
                    'is_binary': self._association_pack.is_binary,
                },
                outputs=['results_file']
            )

        launcher.submit_and_monitor()

        # Collect and process results from subjobs
        completed_models = []
        for result in launcher:
            result_file = InputFileHandler(result['results_file']).get_file_handle()
            with result_file.open('rb') as result_reader:
                gene_dicts = pickle.load(result_reader)
                completed_models.extend(gene_dicts)

        return completed_models


@dxpy.entry_point('run_glm_gene_subjob')
def run_glm_gene_subjob(null_model_dxfile: Dict[str, Any], genotype_dict_dxfile: Dict[str, Any], tarball_prefix: str,
                        is_binary: bool) -> Dict[str, Any]:
    """A subjob to run GLM models for a single tarball prefix.

    This function is executed as a DNAnexus entry point. It deserializes the null model
    and genetic data for a specific tarball prefix, runs gene-level GLMs for all genes
    within that prefix, and then serializes and exports the results.

    :param null_model_dxfile: DNAnexus file ID for the pickled null model.
    :param genotype_dict_dxfile: DNAnexus file ID for the pickled genotype dictionary (pandas DataFrame).
    :param tarball_prefix: The prefix for the current tarball (used as mask_name).
    :param is_binary: Boolean indicating if the phenotype is binary.
    :return: A dictionary containing the DNAnexus file ID for the pickled list of gene results.
    """

    # Load the null model
    null_model_file = InputFileHandler(null_model_dxfile).get_file_handle()
    with null_model_file.open('rb') as null_reader:
        null_model = pickle.load(null_reader)

    # Load the genotype dictionary
    genotype_dict_file = InputFileHandler(genotype_dict_dxfile).get_file_handle()
    genotype_dict = pd.read_pickle(genotype_dict_file)

    # Run linear models for all genes in this chunk
    runner_thread_utility = ThreadUtility()

    for gene in genotype_dict.index.levels[0]:
        runner_thread_utility.launch_job(run_linear_model,
                                         inputs={'linear_model_pack': null_model,
                                                 'genotype_table': genotype_dict,
                                                 'gene': gene,
                                                 'mask_name': tarball_prefix,
                                                 'is_binary': is_binary},
                                         outputs=['gene_dict'])

    runner_thread_utility.submit_and_monitor()

    completed_models = [result['gene_dict'] for result in runner_thread_utility]

    # Save results to a file and export
    results_path = Path("glm_results.pkl")
    with results_path.open('wb') as results_writer:
        pickle.dump(completed_models, results_writer)

    exporter = ExportFileHandler(delete_on_upload=True)
    results_file = exporter.export_files(results_path)

    return {'results_file': results_file}
