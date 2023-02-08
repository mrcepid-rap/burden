from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import dxpy

from general_utilities.association_resources import BGENInformation
from runassociationtesting.association_pack import AssociationPack, ProgramArgs


@dataclass
class BurdenProgramArgs(ProgramArgs):
    association_tarballs: dxpy.DXFile
    tool: str
    run_marker_tests: bool
    bgen_index: dxpy.DXFile
    array_bed_file: dxpy.DXFile
    array_fam_file: dxpy.DXFile
    array_bim_file: dxpy.DXFile
    low_MAC_list: dxpy.DXFile
    sparse_grm: dxpy.DXFile
    sparse_grm_sample: dxpy.DXFile
    bolt_non_infinite: bool
    regenie_smaller_snps: Optional[dxpy.DXFile]


class BurdenAssociationPack(AssociationPack):

    def __init__(self, association_pack: AssociationPack, tarball_prefixes: List[str],
                 bgen_dict: Dict[str, BGENInformation], run_marker_tests: bool, is_bolt_non_infinite: bool,
                 regenie_snps_file: Optional[Path]):

        super().__init__(association_pack.pheno_files, association_pack.inclusion_found,
                         association_pack.exclusion_found, association_pack.additional_covariates_found,
                         association_pack.is_binary, association_pack.sex, association_pack.threads,
                         association_pack.pheno_names,
                         association_pack.found_quantitative_covariates, association_pack.found_categorical_covariates)

        self.tarball_prefixes = tarball_prefixes
        self.bgen_dict = bgen_dict
        self.run_marker_tests = run_marker_tests
        self.is_bolt_non_infinite = is_bolt_non_infinite
        self.regenie_snps_file = regenie_snps_file
