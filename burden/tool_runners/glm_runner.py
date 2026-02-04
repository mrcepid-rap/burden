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
        2. Launching parallel subjobs to run gene-level GLMs for each chromosome.
        3. Processing and annotating the results.
        4. Generating Manhattan plots.
        """

        # 1. Do setup for the linear models & get the null model
        self._logger.info("Loading data and running null Linear Model")
        null_model = linear_model_null(self._association_pack.final_covariates,
                                       self._association_pack.pheno_names[0],
                                       self._association_pack.is_binary,
                                       self._association_pack.found_quantitative_covariates,
                                       self._association_pack.found_categorical_covariates)

        # 2. Launch a subjob for each chromosome
        self._logger.info("Submitting Linear Models to subjobs")
        all_completed_models = self._launch_glm_jobs(null_model)

        # 3. Annotate unformatted results and print final tabular outputs
        self._logger.info("Processing Linear Model results")
        output_path = Path(f'{self._output_prefix}.genes.GLM.stats.tsv')
        self._outputs.extend(process_model_outputs(all_completed_models,
                                                   output_path,
                                                   self._association_pack.tarball_type,
                                                   self._transcripts_table))

        # 4. Generate Manhattan Plots
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

    def _launch_glm_jobs(self, null_model: object) -> List[Dict[str, Any]]:
        """Launch a subjob for each chromosome to run GLM associations.

        :param null_model: A pre-computed null model to use for association testing.
        :return: A list of dictionaries, where each dictionary is a completed model.
        """
        launcher = joblauncher_factory(download_on_complete=True)
        exporter = ExportFileHandler(delete_on_upload=False)

        # First, serialize the null model
        null_model_path = Path("null_model.pkl")
        with null_model_path.open('wb') as null_writer:
            pickle.dump(null_model, null_writer)
        dnanexus_null_model = exporter.export_files(null_model_path)

        # And the transcripts table
        transcripts_table_path = Path("transcripts_table.tsv")
        self._transcripts_table.to_csv(transcripts_table_path, sep='\t', index=True)
        dnanexus_transcripts_table = exporter.export_files(transcripts_table_path)

        # And the list of tarballs
        tarball_prefixes = [p.name for p in self._association_pack.tarball_prefixes]

        for chromosome in self._association_pack.bgen_dict:
            launcher.launch_job(
                function=run_glm_chromosome_subjob,
                inputs={
                    'chromosome': chromosome,
                    'null_model_dxfile': dnanexus_null_model,
                    'transcripts_table_dxfile': dnanexus_transcripts_table,
                    "tarball_prefixes": [str(p) for p in self._association_pack.tarball_prefixes],
                    'tarball_type': str(self._association_pack.tarball_type),
                    'is_binary': self._association_pack.is_binary,
                    'threads': self._association_pack.threads
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


@dxpy.entry_point('run_glm_chromosome_subjob')
def run_glm_chromosome_subjob(chromosome: str, null_model_dxfile: Dict[str, Any],
                              transcripts_table_dxfile: Dict[str, Any], tarball_prefixes: List[str], tarball_type: str,
                              is_binary: bool, threads: int) -> Dict[str, Any]:
    """A subjob to run GLM models for a single chromosome across all tarballs.

    :param chromosome: Chromosome to process.
    :param null_model_dxfile: DNAnexus file ID for the pickled null model.
    :param transcripts_table_dxfile: DNAnexus file ID for the transcripts table.
    :param tarball_prefixes: A list of tarball prefixes to process.
    :param tarball_type: The type of tarballs to load.
    :param is_binary: Boolean indicating if the phenotype is binary.
    :param threads: Number of threads to use for parallel processing.
    :return: A dictionary containing the DNAnexus file ID for the pickled list of gene results.
    """

    # Load the null model
    null_model_file = InputFileHandler(null_model_dxfile).get_file_handle()
    with null_model_file.open('rb') as null_reader:
        null_model = pickle.load(null_reader)

    # Load transcripts table
    transcripts_table_file = InputFileHandler(transcripts_table_dxfile).get_file_handle()
    transcripts_table = pd.read_csv(transcripts_table_file, sep='\t', index_col=0)

    # Get genes for the current chromosome
    genes_on_chrom = transcripts_table[transcripts_table['chrom'] == chromosome].index

    # Load genotype data in parallel for all tarballs
    loader_thread_utility = ThreadUtility(threads=threads, thread_factor=2, incrementor=10)
    for tarball_prefix in tarball_prefixes:
        loader_thread_utility.launch_job(load_linear_model_genetic_data,
                                         inputs={'tarball_prefix': Path(tarball_prefix),
                                                 'tarball_type': tarball_type,
                                                 'bgen_prefix': chromosome},
                                         outputs=['tarball_prefix', 'genotype_dict'])
    loader_thread_utility.submit_and_monitor()

    # Run linear models for all genes in this chunk
    runner_thread_utility = ThreadUtility(threads=threads)
    for result in loader_thread_utility:
        tarball_prefix = result['tarball_prefix']
        genotype_dict = result['genotype_dict']

        genes_to_run = genotype_dict.index.levels[0]

        for gene in genes_to_run:
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
    results_path = Path(f"glm_results.{chromosome}.pkl")
    with results_path.open('wb') as results_writer:
        pickle.dump(completed_models, results_writer)

    exporter = ExportFileHandler(delete_on_upload=True)
    results_file = exporter.export_files(results_path)

    return {'results_file': results_file}
