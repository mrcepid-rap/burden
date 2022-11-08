from typing import Type

import dxpy
from burden import burden_ingester
from burden.burden_association_pack import BurdenProgramArgs, BurdenAssociationPack
from runassociationtesting.module_loader import ModuleLoader
from burden.tool_runners.bolt_runner import BOLTRunner
from burden.tool_runners.glm_runner import GLMRunner
from burden.tool_runners.regenie_runner import REGENIERunner
from burden.tool_runners.saige_runner import SAIGERunner
from burden.tool_runners.staar_runner import STAARRunner
from burden.tool_runners.tool_runner import ToolRunner


class LoadModule(ModuleLoader):

    def __init__(self, output_prefix: str, input_args: str):

        super().__init__(output_prefix, input_args)

    def start_module(self) -> None:

        # Decide which tool we need to run
        current_class = self.check_tools(self.parsed_options.tool)

        # Run the tool – this line just makes an object out of the selected tool we got back from 'check_tools()'. Since
        # every possible tool is a subclass of 'ToolRunner' with a required method of 'run_tool' we should be OK.
        current_tool = current_class(self.association_pack,
                                     self.output_prefix)
        current_tool.run_tool()

        # Retrieve outputs – all tools _should_ append to the outputs object so they can be retrieved here.
        self.set_outputs(current_tool.get_outputs())

    def _load_module_options(self) -> None:

        example_dxfile = 'file-123...'

        self._parser.add_argument('--association_tarballs',
                                  help="Path or hash to list file / single tarball of masks from "
                                       "'mergecollapsevariants'",
                                  type=self.dxfile_input, dest='association_tarballs', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--tool',
                                  help="Select the burden test tool to run (bolt, staar, saige, glm, regenie). "
                                       "Case *MUST* match.",
                                  type=str, dest='tool', required=True,
                                  choices=["bolt", "saige", "staar", "glm", "regenie"])
        self._parser.add_argument('--run_marker_tests',
                                  help="Run per-marker tests for requested tool(s) [true]? Setting to false could "
                                       "DRASTICALLY reduce run-time (particularly if --tool saige). Only changes "
                                       "output for burden tests where tool = saige, regenie, or bolt.",
                                  dest='run_marker_tests', action='store_true')
        self._parser.add_argument('--bgen_index',
                                  help="list of bgen files and associated index/sample/annotation",
                                  type=self.dxfile_input, dest='bgen_index', required=False,
                                  metavar=example_dxfile, default='None')
        self._parser.add_argument('--dosage_index',
                                  help="list of dosage files and associated sample/annotation",
                                  type=self.dxfile_input, dest='dosage_index', required=False,
                                  metavar=example_dxfile, default='None')
        self._parser.add_argument('--array_bed_file',
                                  help="A plink .bed file of genetic data. This file should have been created by "
                                       "mrcepid-buildgrms",
                                  type=self.dxfile_input, dest='array_bed_file', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--array_fam_file',
                                  help="A plink .fam file of genetic data. This file should have been created by "
                                       "mrcepid-buildgrms",
                                  type=self.dxfile_input, dest='array_fam_file', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--array_bim_file',
                                  help="A plink .bim file of genetic data. This file should have been created by "
                                       "mrcepid-buildgrms",
                                  type=self.dxfile_input, dest='array_bim_file', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--low_MAC_list',
                                  help="List of low MAC (max-mac 100) markers in 'bed_file'. This file should have "
                                       "been created by mrcepid-buildgrms",
                                  type=self.dxfile_input, dest='low_MAC_list', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--sparse_grm',
                                  help="A sparse GRM generated by SAIGE step0",
                                  type=self.dxfile_input, dest='sparse_grm', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--sparse_grm_sample',
                                  help="Sparse GRM sample file with row/column names for 'sparse_grm'",
                                  type=self.dxfile_input, dest='sparse_grm_sample', required=True,
                                  metavar=example_dxfile)
        self._parser.add_argument('--bolt_non_infinite',
                                  help="Run bolt in non-infinitesimal mode (bolt flag '--lmmForceNonInf'). WARNING – "
                                       "this can drastically increase runtimes, especially when analysing the "
                                       "entire UKBB dataset.",
                                  dest='bolt_non_infinite', action='store_true')
        self._parser.add_argument('--regenie_smaller_snps',
                                  help="Run REGENIE with the smaller subset of SNPs used for relatedness calculations "
                                       "in Bycroft et al. Input for this parameter is the file ID of the SNP qc file "
                                       "provided by the UKBB "
                                       "[typically located at /Bulk/Genotype Results/Genotype calls/ukb_snp_qc.txt].",
                                  type=self.dxfile_input, dest='regenie_smaller_snps', required=False,
                                  default='None')

    def _parse_options(self) -> BurdenProgramArgs:
        return BurdenProgramArgs(**vars(self._parser.parse_args(self._input_args.split())))

    def _ingest_data(self, parsed_options: BurdenProgramArgs) -> BurdenAssociationPack:
        ingested_data = burden_ingester.BurdenIngestData(parsed_options)
        return ingested_data.get_association_pack()

    # Just defines possible tools usable by this module
    @staticmethod
    def check_tools(input_tool) -> Type[ToolRunner]:

        module_tools = {'bolt': BOLTRunner,
                        'saige': SAIGERunner,
                        'regenie': REGENIERunner,
                        'staar': STAARRunner,
                        'glm': GLMRunner}
        if input_tool in module_tools:
            return module_tools[input_tool]
        else:
            raise dxpy.AppError(f'Tool – {input_tool} – not support. Please try a different input tool!')
