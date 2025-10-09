import csv
from pathlib import Path
from typing import Optional

import dxpy
from general_utilities.import_utils.genetics_loader import GeneticsLoader
from general_utilities.import_utils.import_lib import ingest_wes_bgen, ingest_tarballs, TarballType
from general_utilities.import_utils.module_loader.ingest_data import IngestData, InputFileHandler
from general_utilities.import_utils.import_lib import process_regenie_step_one
from burden.burden_association_pack import BurdenAssociationPack, BurdenProgramArgs


class BurdenIngestData(IngestData):

    def __init__(self, parsed_options: BurdenProgramArgs):
        super().__init__(parsed_options)

        has_regenie_step_one = process_regenie_step_one(parsed_options.regenie_run)

        # Put additional options/covariate processing required by this specific package here
        if len(self.get_association_pack().pheno_names) > 1:
            raise dxpy.AppError('The burden module currently only allows for running one phenotype at a time!')

        tarball_type, tarball_prefixes = ingest_tarballs(parsed_options.association_tarballs)
        if tarball_type is not tarball_type.GENOMEWIDE:
            raise dxpy.AppError('The burden module is not compatible with SNP or GENE masks!')

        # Ingest WES filtered and annotated bgen
        bgen_dict = ingest_wes_bgen(parsed_options.bgen_index)


        loaded_genetics = GeneticsLoader(parsed_options.array_bed_file,
                                         parsed_options.array_fam_file,
                                         parsed_options.array_bim_file,
                                         sample_files=[],
                                         cmd_executor=self._association_pack.cmd_executor,
                                         low_mac_list=parsed_options.low_MAC_list,
                                         sparse_grm=parsed_options.sparse_grm,
                                         sparse_grm_sample=parsed_options.sparse_grm_sample)


        # Put additional covariate processing specific to this module here
        self.set_association_pack(BurdenAssociationPack(self.get_association_pack(),
                                                        tarball_prefixes=tarball_prefixes,
                                                        bgen_dict=bgen_dict,
                                                        run_marker_tests=parsed_options.run_marker_tests,
                                                        is_bolt_non_infinite=parsed_options.bolt_non_infinite,
                                                        has_regenie_step_one=has_regenie_step_one,
                                                        low_mac_list=loaded_genetics.get_low_mac_list(),
                                                        sparse_grm=loaded_genetics.get_sparsematrix(),
                                                        sparse_grm_sample=loaded_genetics.get_sparsematrix_sample(),
                                                        genetic_filename=loaded_genetics.get_filtered_genetic_filename(),
                                                        tarball_type=tarball_type))

