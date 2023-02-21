from burden.tool_runners.tool_runner import ToolRunner
from general_utilities.association_resources import *
from general_utilities.linear_model import linear_model
from general_utilities.linear_model.linear_model import LinearModelResult
from general_utilities.linear_model.proccess_model_output import process_linear_model_outputs
from general_utilities.job_management.thread_utility import *


class GLMRunner(ToolRunner):

    def run_tool(self) -> None:

        # 1. Do setup for the linear models.
        # This will load all variants, genes, and phenotypes into memory to allow for parallelization
        # This function returns a class of type LinearModelPack containing info for running GLMs
        print("Loading data and running null Linear Model")
        null_model = linear_model.linear_model_null(self._association_pack.pheno_names[0],
                                                    self._association_pack.is_binary,
                                                    self._association_pack.found_quantitative_covariates,
                                                    self._association_pack.found_categorical_covariates)

        # 2. Load the tarballs INTO separate genotypes dictionaries
        print("Loading Linear Model genotypes")
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
        print("Submitting Linear Models to threads")
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

        lm_stats_file = open(self._output_prefix + '.lm_stats.tmp', 'w')
        lm_stats_writer = csv.DictWriter(lm_stats_file,
                                         delimiter="\t",
                                         fieldnames=fieldnames,
                                         extrasaction='ignore')

        lm_stats_writer.writeheader()
        for result in thread_utility:
            finished_gene: LinearModelResult = result
            lm_stats_writer.writerow(finished_gene.todict())
        lm_stats_file.close()

        # 5. Annotate unformatted results and print final outputs
        print("Annotating Linear Model results")
        self._outputs.extend(process_linear_model_outputs(self._output_prefix))
