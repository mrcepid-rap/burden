from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import dxpy
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.import_utils.import_lib import BGENInformation
from general_utilities.import_utils.module_loader.association_pack import AssociationPack, ProgramArgs


@dataclass
class BurdenProgramArgs(ProgramArgs):
    association_tarballs: InputFileHandler
    tool: str
    run_marker_tests: bool
    bgen_index: InputFileHandler
    array_bed_file: InputFileHandler
    array_fam_file: InputFileHandler
    array_bim_file: InputFileHandler
    low_MAC_list: InputFileHandler
    sparse_grm: InputFileHandler
    sparse_grm_sample: InputFileHandler
    bolt_non_infinite: bool
    regenie_smaller_snps: Optional[InputFileHandler]

    def __post_init__(self):
        """@dataclass automatically calls this method after calling its own __init__().

        This is required in the subclass because dataclasses do not call the __init__ of their super o.0

        """

        self._check_opts()

    def _check_opts(self):
        pass


class BurdenAssociationPack(AssociationPack):

    def __init__(self, association_pack: AssociationPack, tarball_prefixes: List[str],
                 bgen_dict: Dict[str, BGENInformation], run_marker_tests: bool, is_bolt_non_infinite: bool,
                 regenie_snps_file: Optional[Path], tarball_chunks: List[str], low_mac_list: Path, sparse_grm: Path,
                 sparse_grm_sample: Path, genetic_filename: str):
        super().__init__(association_pack.is_binary, association_pack.sex, association_pack.threads,
                         association_pack.pheno_names, association_pack.ignore_base_covariates,
                         association_pack.found_quantitative_covariates, association_pack.found_categorical_covariates,
                         association_pack.cmd_executor, association_pack.final_covariates, association_pack.inclusion_samples,
                         association_pack.exclusion_samples)

        self.tarball_prefixes = tarball_prefixes
        self.tarball_chunks = tarball_chunks
        self.bgen_dict = bgen_dict
        self.run_marker_tests = run_marker_tests
        self.is_bolt_non_infinite = is_bolt_non_infinite
        self.regenie_snps_file = regenie_snps_file
        self.low_mac_list = low_mac_list
        self.sparse_grm = sparse_grm
        self.sparse_grm_sample = sparse_grm_sample
        self.genetic_filename = genetic_filename
