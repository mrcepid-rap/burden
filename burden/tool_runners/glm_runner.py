"""Execute gene-level GLM runs and downstream plotting for burden analyses."""
import os
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
    load_linear_model_genetic_data, run_linear_model, TarballType
from general_utilities.linear_model.proccess_model_output import process_model_outputs
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter

from burden.tool_runners.tool_runner import ToolRunner


class GLMRunner(ToolRunner):
    """Run the GLM workflow over all tarball prefixes and available genes."""

    def run_tool(self) -> None:
        """Load input data, fit models, persist results, and generate plots."""

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

        # Filter out results for genes that are not in the transcripts table
        valid_genes = set(self._transcripts_table.index)
        all_completed_models = [m for m in all_completed_models if m.ENST in valid_genes]

        if all_completed_models:
            output_path = Path(f'{self._output_prefix}.genes.GLM.stats.tsv')
            self._outputs.extend(process_model_outputs(all_completed_models,
                                                       output_path,
                                                       self._association_pack.tarball_type,
                                                       self._transcripts_table))

            # 4. Generate Manhattan Plots
            plot_dir = Path(f'{self._output_prefix}_plots')
            plot_dir.mkdir(parents=True, exist_ok=True)
            self._outputs.append(plot_dir)

            stats_gz = Path(f'{self._output_prefix}.genes.GLM.stats.tsv.gz')
            if stats_gz.exists():
                glm_table = pd.read_csv(stats_gz, sep='\t')

                if not glm_table.empty:
                    mask_values = glm_table['MASK'].value_counts().index
                    maf_values = glm_table['MAF'].value_counts().index

                    for mask in mask_values:
                        for maf in maf_values:
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

                            plots = manhattan_plotter.plot()
                            if plots:
                                plots[0].rename(plot_dir / f'{mask}.{maf}.genes.GLM.png')
        else:
            self._logger.warning("No valid models found after filtering.")

    def _launch_glm_jobs(self, null_model: object) -> List[Dict[str, Any]]:
        """Launch a subjob for each chromosome to run GLM associations."""
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

        for chromosome in self._association_pack.bgen_dict:

            # tarball BGEN files
            bgen_filename = [f for f in Path('.').glob(f'*{chromosome}.BOLT.bgen') if not f.name.startswith('._')]
            bgen_index_filename = [f for f in Path('.').glob(f'*{chromosome}.BOLT.bgen.bgi') if
                                   not f.name.startswith('._')]
            bgen_sample_filename = [f for f in Path('.').glob(f'*{chromosome}.BOLT.sample') if
                                    not f.name.startswith('._')]

            # Only launch job if BOLT bgen file exists
            if not bgen_filename:
                self._logger.info(f"No BOLT bgen file found for chromosome {chromosome}, skipping job launch.")
                continue

            # export BOLT bgen files
            bolt_bgen = exporter.export_files(bgen_filename)
            bolt_bgen_index = exporter.export_files(bgen_index_filename)
            bolt_bgen_sample = exporter.export_files(bgen_sample_filename)

            launcher.launch_job(
                function=run_glm_chromosome_subjob,
                inputs={
                    'chromosome': chromosome,
                    'null_model_dxfile': dnanexus_null_model,
                    'transcripts_table_dxfile': dnanexus_transcripts_table,
                    "tarball_prefixes": [str(p) for p in self._association_pack.tarball_prefixes],
                    'tarball_type': self._association_pack.tarball_type.value,
                    'is_binary': self._association_pack.is_binary,
                    'threads': self._association_pack.threads,
                    'bolt_bgen': bolt_bgen,
                    'bolt_bgen_index': bolt_bgen_index,
                    'bolt_bgen_sample': bolt_bgen_sample
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
                              bolt_bgen: List[Any], bolt_bgen_index: List[Any], bolt_bgen_sample: List[Any],
                              is_binary: bool, threads: int) -> Dict[str, Any]:

    # Download BOLT Data (Fix: iterate because input is a list)
    for f in bolt_bgen:
        InputFileHandler(f, download_now=True).get_file_handle()
    for f in bolt_bgen_index:
        InputFileHandler(f, download_now=True).get_file_handle()
    for f in bolt_bgen_sample:
        InputFileHandler(f, download_now=True).get_file_handle()

    # Load the null model
    null_model_file = InputFileHandler(null_model_dxfile).get_file_handle()
    with null_model_file.open('rb') as null_reader:
        null_model = pickle.load(null_reader)

    # Load transcripts table
    transcripts_table_file = InputFileHandler(transcripts_table_dxfile).get_file_handle()
    transcripts_table = pd.read_csv(transcripts_table_file, sep='\t', index_col=0)

    # Convert the string back to the Enum object required by load_linear_model_genetic_data
    tarball_type_enum = TarballType(tarball_type)

    # Load genotype data in parallel for all tarballs
    loader_thread_utility = ThreadUtility(threads=threads, thread_factor=2, incrementor=10)
    for tarball_prefix in tarball_prefixes:
        loader_thread_utility.launch_job(load_linear_model_genetic_data,
                                         inputs={'tarball_prefix': Path(tarball_prefix),
                                                 'tarball_type': tarball_type_enum,
                                                 'bgen_prefix': chromosome},
                                         outputs=['tarball_prefix', 'genotype_dict'])
    loader_thread_utility.submit_and_monitor()

    # Run linear models for all genes in this chunk
    runner_thread_utility = ThreadUtility(threads=threads)
    for result in loader_thread_utility:
        current_tarball_prefix = result['tarball_prefix']
        genotype_dict = result['genotype_dict']

        genes_to_run = genotype_dict.index.levels[0]

        for gene in genes_to_run:
            runner_thread_utility.launch_job(run_linear_model,
                                             inputs={'linear_model_pack': null_model,
                                                     'genotype_table': genotype_dict,
                                                     'gene': gene,
                                                     'mask_name': current_tarball_prefix,
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