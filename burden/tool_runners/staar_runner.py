import json
from pathlib import Path
from typing import Any, List, Dict

import dxpy
import pandas as pd
from general_utilities.bgen_utilities.genotype_matrix import generate_csr_matrix_from_bgen
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import LOGGER, TarballType
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.joblauncher_factory import joblauncher_factory
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.linear_model.proccess_model_output import process_model_outputs
from general_utilities.linear_model.staar_model import load_staar_genetic_data, staar_genes, staar_null, \
    STAARModelResult
from general_utilities.plot_lib.manhattan_plotter import ManhattanPlotter
from scipy.io import mmwrite
from dataclasses import asdict
from burden.tool_runners.tool_runner import ToolRunner
import gc


class STAARRunner(ToolRunner):
    """Coordinate STAAR null models, per-gene burden tests, and downstream plots."""

    def run_tool(self) -> None:
        """Execute the STAAR pipeline across phenotypes, tarball masks, and chromosomes."""

        # 1. Run the STAAR NULL model
        self._logger.info("Running STAAR Null Model...")
        if not Path('phenotype.STAAR_null.rds').exists():
            # 1. Run the STAAR NULL model
            self._logger.info("Running STAAR Null Model(s)...")
            thread_utility = ThreadUtility(self._association_pack.threads,
                                           thread_factor=1)

            first_chrom = next(iter(self._association_pack.bgen_dict))
            sample_path = self._association_pack.bgen_dict[first_chrom]["sample"].get_file_handle()

            # Ensure both columns are strings for merging
            sample = pd.read_csv(sample_path, sep=r"\s+", header=0, dtype={'ID_2': str})
            sample = sample.drop(columns=["sex"], errors="ignore")
            sample = sample.iloc[1:].reset_index(drop=True)

            covar = pd.read_csv(self._association_pack.final_covariates, sep=' ', header=0, dtype={'IID': str})

            merged = sample.merge(
                covar,
                how="left",
                left_on="ID_2",
                right_on="IID"
            )

            # 6. Drop only these (sex already removed)
            merged = merged.drop(columns=["ID_1", "missing"], errors="ignore")

            # 7. Sort by ID_2
            merged = merged.sort_values("ID_2")

            pd.set_option('display.max_columns', None)

            print("=== MERGED HEAD ===")
            print(merged.head())

            # 8. Write merged covariates to file for STAAR Null
            merged_cov_path = Path("merged_covariates_for_staar.tsv")

            # Remove rows with missing phenotype or any covariates that will be used in the model
            # This ensures the sample list matches what R will actually use
            for phenoname in self._association_pack.pheno_names:
                if phenoname in merged.columns:
                    merged = merged.dropna(subset=[phenoname])

            # Also drop samples with missing covariates that will be used
            required_cols = ['age', 'age_squared', 'batch']
            if self._association_pack.sex == 2:
                required_cols.append('sex')
            for i in range(1, 11):
                required_cols.append(f'PC{i}')
            required_cols.extend(self._association_pack.found_quantitative_covariates)
            required_cols.extend(self._association_pack.found_categorical_covariates)

            # Only keep columns that exist
            required_cols = [col for col in required_cols if col in merged.columns]
            merged = merged.dropna(subset=required_cols)

            self._logger.info(f"After removing samples with missing data: {len(merged)} samples remain")

            merged.to_csv(merged_cov_path, sep="\t", index=False)


            for phenoname in self._association_pack.pheno_names:
                thread_utility.launch_job(function=staar_null,
                                          inputs={
                                              'phenofile': str(merged_cov_path),
                                              'phenotype': phenoname,
                                              'is_binary': self._association_pack.is_binary,
                                              'ignore_base': self._association_pack.ignore_base_covariates,
                                              'found_quantitative_covariates': self._association_pack.found_quantitative_covariates,
                                              'found_categorical_covariates': self._association_pack.found_categorical_covariates,
                                              'sex': self._association_pack.sex,
                                              'sparse_kinship_file': self._association_pack.sparse_grm,
                                              'sparse_kinship_samples': self._association_pack.sparse_grm_sample
                                          },
                                          outputs=['staar_null_model']
                                          )
            thread_utility.submit_and_monitor()
        else:
            self._logger.info("STAAR null model already exists, skipping...")

        # Filter STAAR samples tables to match the samples in the null model
        self._logger.info("Filtering STAAR samples tables to match null model samples...")

        # Read the samples that were used in the null model
        merged_cov_path = Path("merged_covariates_for_staar.tsv")
        if merged_cov_path.exists():
            merged_samples = pd.read_csv(merged_cov_path, sep='\t')
            null_model_samples = set(merged_samples['ID_2'].astype(str))

            # Now filter each STAAR samples table
            for tarball_prefix in self._association_pack.tarball_prefixes:
                for chromosome in self._association_pack.bgen_dict:
                    samples_path = Path(f"{tarball_prefix}.{chromosome}.STAAR.samples_table.tsv")

                    if samples_path.exists():
                        # Read the original STAAR samples table
                        staar_samples_df = pd.read_csv(samples_path, sep='\t')
                        original_count = len(staar_samples_df)

                        # Filter to only keep samples in the null model
                        staar_samples_df = staar_samples_df[
                            staar_samples_df['sampID'].astype(str).isin(null_model_samples)
                        ]
                        filtered_count = len(staar_samples_df)

                        # Overwrite the STAAR samples table with the filtered version
                        staar_samples_df.to_csv(samples_path, sep='\t', index=False)

                        self._logger.info(
                            f"Filtered {samples_path.name}: {original_count} → {filtered_count} samples"
                        )

        # 2. Run the actual per-gene association tests
        self._logger.info("Running STAAR masks * chromosomes...")

        # set the launcher
        launcher = joblauncher_factory(download_on_complete=True)

        # set the exporter
        exporter = ExportFileHandler(delete_on_upload=False)

        for phenoname in self._association_pack.pheno_names:
            for tarball_prefix in self._association_pack.tarball_prefixes:
                for chromosome in self._association_pack.bgen_dict:

                    # check that the STAAR files exist for this tarball / chromosome
                    samples_path = Path(f"{tarball_prefix}.{chromosome}.STAAR.samples_table.tsv")
                    variants_path = Path(f"{tarball_prefix}.{chromosome}.STAAR.variants_table.tsv")
                    # Skip if this chromosome was not produced for this tarball
                    if not samples_path.exists() or not variants_path.exists():
                        self._logger.info(f"Skipping {tarball_prefix} / {chromosome}: no STAAR files present")
                        continue

                    # get the STAAR genetic data
                    staar_data = load_staar_genetic_data(
                        tarball_prefix=str(tarball_prefix),
                        bgen_prefix=chromosome
                    )
                    # send the genetic information to each subjob
                    subset_staar_data = {chromosome: staar_data[chromosome]}
                    chunk_json_path = Path(f"{tarball_prefix}.{chromosome}.staar_chunk.json")
                    # Serialize only the current chromosome so each worker downloads a minimal payload.
                    with chunk_json_path.open("w") as f:
                        json.dump(subset_staar_data, f, default=lambda o: list(o) if isinstance(o, set) else o)

                    # export the files we need
                    null_model = exporter.export_files(f'{phenoname}.STAAR_null.rds')
                    staar_samples = exporter.export_files(f'{tarball_prefix}.{chromosome}.STAAR.samples_table.tsv')
                    variants_table = exporter.export_files(
                        f'{tarball_prefix}.{chromosome}.STAAR.variants_table.tsv')
                    chunk_file = exporter.export_files(chunk_json_path)
                    transcripts_table = Path("transcripts_table.tsv")
                    self._transcripts_table.to_csv(transcripts_table, sep='\t', index=True)
                    transcripts_table = exporter.export_files(transcripts_table)

                    # set the chunk that we are working with
                    working_chunk = self._association_pack.bgen_dict[chromosome]

                    launcher.launch_job(
                        function=multithread_staar_burden,
                        inputs={
                            'tarball_prefix': tarball_prefix,
                            'chromosome': chromosome,
                            'phenoname': phenoname,
                            'staar_null_model': null_model,
                            'bgen_file': working_chunk['bgen'].get_input_str(),
                            'bgen_index': working_chunk['index'].get_input_str(),
                            'bgen_sample': working_chunk['sample'].get_input_str(),
                            'variants_table': variants_table,
                            'staar_samples': staar_samples,
                            'chunk_file': chunk_file,
                            'threads': self._association_pack.threads
                        },
                        outputs=['output_model']
                    )
        launcher.submit_and_monitor()

        # Gather all the results
        completed_staar_chunks = []
        for result in launcher:
            # Annotate STAAR outputs
            output_path = Path(f'{self._output_prefix}.genes.STAAR.stats.tsv')
            self._outputs.extend(process_model_outputs(input_models=result['output_model'],
                                  output_path=output_path,
                                  tarball_type=self._association_pack.tarball_type,
                                  transcripts_table=self._transcripts_table))

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
                                                     p_column='p_val_burden',
                                                     csq_column='MASK',
                                                     maf_column='cMAC', gene_symbol_column='SYMBOL',
                                                     clumping_distance=1,
                                                     maf_cutoff=30,
                                                     sig_threshold=1E-6)

                manhattan_plotter.plot()[0].rename(plot_dir / f'{mask}.{maf}.genes.STAAR.png')

@dxpy.entry_point('multithread_staar_burden')
def multithread_staar_burden(tarball_prefix: str, chromosome: str, phenoname: str, staar_null_model: dict,
                             bgen_file: str,
                             bgen_index: str, bgen_sample: str, variants_table: dict, staar_samples: dict,
                             chunk_file: dict, threads: int) -> Dict[str, List[Dict[str, Any]]]:
    """
    Run the STAAR gene tests for a single tarball/chromosome chunk inside the worker environment.

    :param tarball_prefix: Mask identifier used to label inputs/outputs.
    :param chromosome: Chromosome identifier for the current job.
    :param phenoname: Phenotype currently being evaluated.
    :param staar_null_model: Remote path to the serialized null model.
    :param bgen_file: Remote BGEN file pointer delivered by the job framework.
    :param bgen_index: Remote BGEN index pointer delivered by the job framework.
    :param bgen_sample: Remote BGEN sample file pointer delivered by the job framework.
    :param variants_table: Remote variants metadata table.
    :param staar_samples: Remote samples metadata table.
    :param chunk_file: Remote JSON chunk describing the subset of genes to evaluate.
    :param threads: Number of threads to use inside this worker.

    :return: Dictionary containing list of STAAR results.
    """
    # load our VM environment
    cmd_executor = build_default_command_executor()
    null_model = InputFileHandler(staar_null_model).get_file_handle()
    staar_samples = InputFileHandler(staar_samples).get_file_handle()
    staar_variants = InputFileHandler(variants_table).get_file_handle()
    chunk_file = InputFileHandler(chunk_file).get_file_handle()

    with open(chunk_file, "r") as f:
        staar_data = json.load(f)

    # download genetic data
    bgen_path = InputFileHandler(bgen_file).get_file_handle()
    index_path = InputFileHandler(bgen_index).get_file_handle()
    sample_path = InputFileHandler(bgen_sample).get_file_handle()

    # Limit concurrency per worker so that R-based jobs do not overwhelm the VM.
    thread_utility = ThreadUtility(threads=1)

    for gene, info in staar_data[chromosome].items():
        n_variants = len(info.get("vars", []))
        if n_variants < 2:
            LOGGER.warning(f"Skipping {gene}: fewer than 2 variants ({n_variants})")
            continue

        # generate a csr matrix from the bgen files
        matrix, summary_dict = generate_csr_matrix_from_bgen(
            bgen_path=bgen_path,
            sample_path=sample_path,
            variant_filter_list=staar_data[chromosome][gene]['vars'],
            chromosome=staar_data[chromosome][gene]['chrom'],
            start=staar_data[chromosome][gene]['min'],
            end=staar_data[chromosome][gene]['max'],
            should_collapse_matrix=False
        )

        # export matrix to file
        matrix_file = f"{tarball_prefix}.{chromosome}.{gene}.STAAR.mtx"
        mmwrite(matrix_file, matrix)

        # Clean up matrix from memory immediately
        del matrix
        del summary_dict
        gc.collect()

        thread_utility.launch_job(
            function=staar_genes,
            inputs={
                'staar_null_path': null_model,
                'pheno_name': phenoname,
                'gene': gene,
                'mask_name': tarball_prefix,
                'staar_matrix': matrix_file,
                'staar_samples': staar_samples,
                'staar_variants': staar_variants,
                'out_dir': Path('.'),
            },
            outputs=['staar_result']
        )

    thread_utility.submit_and_monitor()

    # Collect results without keeping them all in memory at once
    completed_staar_files = []
    for result in thread_utility:
        staar_result = result["staar_result"]
        # Convert to dict and append
        completed_staar_files.append(asdict(staar_result))

        # Clean up after processing each result
        del result
        del staar_result
        gc.collect()

    return {'output_model': completed_staar_files}