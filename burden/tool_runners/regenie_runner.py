import re
from pathlib import Path
from typing import List, Any, Tuple, Dict

import dxpy
from general_utilities.association_resources import define_covariate_string, \
    define_field_names_from_tarball_prefix, bgzip_and_tabix
    bgzip_and_tabix, define_field_names_from_tarball_prefix
from general_utilities.import_utils.file_handlers.export_file_handler import ExportFileHandler
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import download_bgen_file, LOGGER
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.joblauncher_factory import joblauncher_factory

# TODO: Implement multi-phenotype testing for REGENIE
class REGENIERunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Run step 1 of regenie
        # if step one has already been run (if the files exist and are not empty)
        if not self._association_pack.has_regenie_step_one:
            self._logger.info("Running REGENIE step 1")
            regenie_step1_log = self._run_regenie_step_one()
            # Add the step1 files to output, so we can use later if need-be:
            self._outputs.extend([Path('fit_out_pred.list'),
                                  Path('fit_out_1.loco'),
                                  regenie_step1_log])
        else:
            self._logger.info("REGENIE step 1 files already exist, skipping...")

        # # 2. Prep bgen files for a run:
        self._logger.info("Downloading and filtering raw bgen files")
        thread_utility = ThreadUtility(self._association_pack.threads,
                                       thread_factor=4)

        for chromosome in self._association_pack.bgen_dict:
        self._logger.info("Preparing to launch subjobs for running REGENIE step 2 on separate VMs...")
            # IDENTICAL to that done for BOLT.
        # set the output
        all_step2_outputs = []

        # --- CHANGED: Multithread by chromosome ---
        thread_utility = ThreadUtility(thread_factor=1)
        # 4. Run step 2 of regenie
            thread_utility.launch_job(
                function=self._multithread_step2,
                inputs=(chromosome,),
                outputs=["step2_outputs"]
            )
        thread_utility.submit_and_monitor()
                if Path(f'{tarball_prefix}.{chromosome}.REGENIE.annotationFile.txt').exists():
        for result in thread_utility:
            all_step2_outputs.extend(result["step2_outputs"])
        # --- END CHANGED ---
                                              chromosome=chromosome)

        # Gather preliminary results from step 2:
        self._logger.info("Gathering REGENIE mask-based results...")
        completed_gene_tables = []
        regenie_step2_genes_log = Path(f'{self._output_prefix}.REGENIE_step2_genes.log')
            for result in all_step2_outputs:
                tarball_prefix = result['tarball_prefix']
                finished_chromosome = result['finished_chromosome']
                current_log = result['current_log']
                completed_gene_tables.append(self._process_regenie_output(tarball_prefix,
                                                                          finished_chromosome))

                # Write a header for each logfile
                regenie_step2_genes_writer.write(f'{tarball_prefix + "-" + finished_chromosome:{"-"}^{50}}')

                with current_log.open('r') as current_log_reader:
                    for line in current_log_reader:
                        regenie_step2_genes_writer.write(line)

        self._outputs.append(regenie_step2_genes_log)
            self._outputs.append(regenie_step2_genes_log)
        # Always add predictions...
            # Always add predictions...
            self._outputs.append(Path('fit_out_pred.list'))
            self._outputs.extend(Path('./').glob('fit_out_*.loco'))
        # 5. Process outputs
            # 5. Process outputs
            self._logger.info("Processing REGENIE outputs...")
            self._outputs.extend(self._annotate_regenie_output(completed_gene_tables=completed_gene_tables))
    def _run_regenie_step_one(self) -> Path:

        # And generate a SNP list for the --extract parameter of REGENIE, while considering SNPs from
        # the regenie_smaller_snps input parameter (if provided). plink2 order of operations:
        # 1. Select variants from --extract (if present)
        # 2. THEN filter based on max/min AC (mac/max-mac)
        max_mac = (self._sample_count * 2) - 100
        cmd = f'plink2 --bfile /test/{self._association_pack.genetic_filename} ' \
        cmd = f"plink2 --bed {self._association_pack.genetic_filename}.bed " \
              f"--bim {self._association_pack.genetic_filename}.bim " \
              f"--fam {self._association_pack.genetic_filename}.fam " \
              f"--min-ac 100 " \

        with open('plink_out.txt', 'r') as plink_out:
              f'--out REGENIE_extract'
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd)
                if found_snp_count is not None:
        with open('REGENIE_extract.snplist', 'r') as plink_out:
            plink_out.close()

        cmd = 'regenie ' \
              '--step 1 ' \
              f'--bed /test/{self._association_pack.genetic_filename} ' \
              '--extract /test/REGENIE_extract.snplist ' \
              f'--covarFile /test/{self._association_pack.final_covariates} ' \
              f'--phenoFile /test/{self._association_pack.final_covariates} ' \
              f'--bed {self._association_pack.genetic_filename} ' \
              '--extract REGENIE_extract.snplist ' \
              f'--covarFile {self._association_pack.final_covariates} ' \
              f'--phenoFile {self._association_pack.final_covariates} ' \
              f'--phenoCol {self._association_pack.pheno_names[0]} '

              '--out fit_out ' \
                                       self._association_pack.found_categorical_covariates,
                                       self._association_pack.is_binary,
                                       add_array=False,
                                       ignore_base=self._association_pack.ignore_base_covariates)

        regenie_log = Path(f'{self._output_prefix}.REGENIE_step1.log')
        self._association_pack.cmd_executor.run_cmd_on_docker(cmd, stdout_file=regenie_log)
        return regenie_log

    def _run_regenie_step_two(self, tarball_prefix: str, chromosome: str) -> Tuple[str, str, Path]:


        cmd = f'regenie ' \
              f'--step 2 ' \
    def _multithread_step2(self, chromosome) -> List[Dict[str, Any]]:
        """
        A function to run REGENIE step 2 on a single chromosome in a multithreaded manner

        :param chromosome: The chromosome/chunk to run
        :return: A list of dictionaries containing the tarball prefix, finished chromosome, and log file
        """

        # set the launcher
        launcher = joblauncher_factory()

        # files to include
        samples_include = Path('SAMPLES_Include.txt')
        fit_out_pred = Path('fit_out_pred.list')
        fit_out_loco = Path('fit_out_1.loco')
              f'--ref-first ' \
        # set the exporter
        exporter = ExportFileHandler()

        # make a list of all annotation files for this chromosome
        anno_files = list(Path('.').glob(f'*.{chromosome}.REGENIE.annotationFile.txt'))
        # make a list of the mask files for this chromosome
        mask_files = list(Path('.').glob(f'*.{chromosome}.REGENIE.maskfile.txt'))
        # make a list of the setlist files for this chromosome
        setlist_files = list(Path('.').glob(f'*.{chromosome}.REGENIE.setListFile.txt'))
        regenie_log = Path(f'{tarball_prefix}.{chromosome}.REGENIE_genes.log')
        # export the files to DX for each subjob
        samples_include = exporter.export_files(samples_include)
        fit_out_pred = exporter.export_files(fit_out_pred)
        fit_out_loco = exporter.export_files(fit_out_loco)
        phenotype_file = exporter.export_files(self._association_pack.final_covariates)
        anno_files = [exporter.export_files(af) for af in anno_files]
        mask_files = [exporter.export_files(mf) for mf in mask_files]
        setlist_files = [exporter.export_files(sf) for sf in setlist_files]
                                    comment='#')
        print(self._association_pack.bgen_dict[chromosome]['bgen'].get_input_str())
        print(self._association_pack.bgen_dict[chromosome]['sample'].get_input_str())
        print(chromosome)
        print(self._association_pack.tarball_prefixes)
        print(samples_include)
        print(self._association_pack.final_covariates)
        print(self._association_pack.pheno_names[0])
        print(fit_out_pred)
        print(fit_out_loco)
        print(anno_files)
        print(mask_files)
        print(setlist_files)
        print(self._association_pack.found_quantitative_covariates)
        print(self._association_pack.found_categorical_covariates)
        print(self._association_pack.is_binary)
        print(self._association_pack.ignore_base_covariates)

        launcher.launch_job(function=run_regenie_step2,
                            inputs={
                                "bgen_file": {'$dnanexus_link': self._association_pack.bgen_dict[chromosome]['bgen'].get_input_str()},
                                "bgen_sample": {'$dnanexus_link': self._association_pack.bgen_dict[chromosome]['sample'].get_input_str()},
                                "bgen_index": {'$dnanexus_link': self._association_pack.bgen_dict[chromosome]['index'].get_input_str()},
                                "chromosome": chromosome,
                                "tarball_prefixes": self._association_pack.tarball_prefixes,
                                "samples_include": samples_include,
                                "covariate_file": phenotype_file,
                                "pheno_file": phenotype_file,
                                "pheno_column": self._association_pack.pheno_names[0],
                                "fit_out_pred": fit_out_pred,
                                "fit_out_loco": fit_out_loco,
                                "annotation_file": anno_files,
                                "mask_file": mask_files,
                                "setlist_file": setlist_files,
                                "found_quantitative_covariates": self._association_pack.found_quantitative_covariates,
                                "found_categorical_covariates": self._association_pack.found_categorical_covariates,
                                "is_binary": self._association_pack.is_binary,
                                "ignore_base_covariates": self._association_pack.ignore_base_covariates,
                            },
                            outputs=[
                                "output"
                            ]
                            )

        launcher.submit_and_monitor()

        step2_outputs = []
        for result in launcher:
            # result["output"] is already a list of dicts
            for r in result["output"]:
                # download subjob outputs to local machine
                InputFileHandler(r["current_log"], download_now=True).get_file_handle()
                InputFileHandler(r["regenie_output"], download_now=True).get_file_handle()
                step2_outputs.append(r)

        return step2_outputs

            # Sort just in case
            regenie_table = regenie_table.sort_values(by=['chrom', 'start', 'end'])

            for mask in regenie_table['MASK'].value_counts().index:

                for maf in regenie_table['MAF'].value_counts().index:
                    # To note on the below: I use SYMBOL for the id_column parameter below because ENST is the
                    # index and I don't currently have a way to pass the index through to the Plotter methods...
                    manhattan_plotter = ManhattanPlotter(self._association_pack.cmd_executor,
                                                         regenie_table.query(f'MASK == "{mask}" & MAF == "{maf}"'),
                                                         chrom_column='chrom', pos_column='start',
                                                         alt_column=None,
                                                         id_column='SYMBOL', p_column='PVALUE',
                                                         csq_column='MASK',
                                                         maf_column='A1FREQ', gene_symbol_column='SYMBOL',
                                                         clumping_distance=1,
                                                         maf_cutoff=30 / (regenie_table['N'].max() * 2),
                                                         sig_threshold=1E-6)

                    manhattan_plotter.plot()[0].rename(plot_dir / f'{mask}.{maf}.genes.REGENIE.png')
            # Write to disk
            regenie_table.to_csv(path_or_buf=gene_out, index=False, sep="\t", na_rep='NA')

        # And bgzip and tabix...
        outputs.extend(bgzip_and_tabix(regenie_gene_out, comment_char='c', sequence_row=1, begin_row=2, end_row=3))

        return outputs

