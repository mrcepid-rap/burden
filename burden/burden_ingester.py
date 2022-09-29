import tarfile

from os.path import exists
from pathlib import Path
from typing import Dict, Optional

from burden.burden_association_pack import BurdenAssociationPack, BGENInformation, \
    BurdenProgramArgs
from runassociationtesting.ingest_data import *


class BurdenIngestData(IngestData):

    def __init__(self, parsed_options: BurdenProgramArgs):
        super().__init__(parsed_options)

        # Put additional options/covariate processing required by this specific package here
        if len(self.get_association_pack().pheno_names) > 1:
            raise dxpy.AppError('The burden module currently only allows for running one phenotype at a time!')

        is_snp_tar, is_gene_tar, tarball_prefixes = self._ingest_tarballs(parsed_options.association_tarballs)
        if is_snp_tar or is_gene_tar:
            raise dxpy.AppError('The burden module is not compatible with SNP or GENE masks!')

        bgen_dict = self._ingest_bgen(parsed_options.bgen_index)

        self._ingest_genetic_data(parsed_options.array_bed_file,
                                  parsed_options.array_fam_file,
                                  parsed_options.array_bim_file,
                                  parsed_options.low_MAC_list,
                                  parsed_options.sparse_grm,
                                  parsed_options.sparse_grm_sample)

        self._generate_filtered_genetic_data()
        regenie_snps_file = self._process_regenie_snps(parsed_options.regenie_smaller_snps)

        # Put additional covariate processing specific to this module here
        self.set_association_pack(BurdenAssociationPack(self.get_association_pack(),
                                                        tarball_prefixes, bgen_dict, parsed_options.run_marker_tests,
                                                        parsed_options.bolt_non_infinite, regenie_snps_file))

    # Need to grab the tarball file for associations...
    # This was generated by the applet mrcepid-collapsevariants
    # Ingest the list file into this AWS instance
    @staticmethod
    def _ingest_tarballs(association_tarballs: dxpy.DXFile) -> Tuple[bool, bool, List[str]]:

        is_snp_tar = False
        is_gene_tar = False
        tarball_prefixes = []
        if '.tar.gz' in association_tarballs.describe()['name']:
            # likely to be a single tarball, download, check, and extract:
            tarball_name = association_tarballs.describe()['name']
            dxpy.download_dxfile(association_tarballs, tarball_name)
            if tarfile.is_tarfile(tarball_name):
                tarball_prefix = tarball_name.replace(".tar.gz", "")
                tarball_prefixes.append(tarball_prefix)
                tar = tarfile.open(tarball_name, "r:gz")
                tar.extractall()
                if exists(tarball_prefix + ".SNP.BOLT.bgen"):
                    is_snp_tar = True
                elif exists(tarball_prefix + ".GENE.BOLT.bgen"):
                    is_gene_tar = True
            else:
                raise dxpy.AppError(f'Provided association tarball ({association_tarballs.describe()["id"]}) '
                                    f'is not a tar.gz file')
        else:
            # Likely to be a list of tarballs, download and extract...
            dxpy.download_dxfile(association_tarballs, "tarball_list.txt")
            with open("tarball_list.txt", "r") as tarball_reader:
                for association_tarball in tarball_reader:
                    association_tarball = association_tarball.rstrip()
                    tarball = dxpy.DXFile(association_tarball)
                    tarball_name = tarball.describe()['name']
                    dxpy.download_dxfile(tarball, tarball_name)

                    # Need to get the prefix on the tarball to access resources within:
                    # All files within SHOULD have the same prefix as this file
                    tarball_prefix = tarball_name.rstrip('.tar.gz')
                    tarball_prefixes.append(tarball_prefix)
                    tar = tarfile.open(tarball_name, "r:gz")
                    tar.extractall()
                    if exists(tarball_prefix + ".SNP.BOLT.bgen"):
                        raise dxpy.AppError(f'Cannot run masks from a SNP list ({association_tarballs.describe()["id"]}) '
                                            f'when running tarballs as batch...')
                    elif exists(tarball_prefix + ".GENE.BOLT.bgen"):
                        raise dxpy.AppError(f'Cannot run masks from a GENE list ({association_tarballs.describe()["id"]}) '
                                            f'when running tarballs as batch...')

        return is_snp_tar, is_gene_tar, tarball_prefixes

    # Grab the entire WES variant data in bgen format
    @staticmethod
    def _ingest_bgen(bgen_index: dxpy.DXFile) -> Dict[str, BGENInformation]:

        # Ingest the INDEX of bgen files:
        dxpy.download_dxfile(bgen_index.get_id(), "bgen_locs.tsv")
        # and load it into a dict:
        os.mkdir("filtered_bgen/")  # For downloading later...
        bgen_index_csv = csv.DictReader(open("bgen_locs.tsv", "r"), delimiter="\t")
        bgen_dict = {}
        for line in bgen_index_csv:
            bgen_dict[line['chrom']] = {'index': line['bgen_index_dxid'],
                                        'sample': line['sample_dxid'],
                                        'bgen': line['bgen_dxid'],
                                        'vep': line['vep_dxid']}
        return bgen_dict

    @staticmethod
    def _ingest_genetic_data(bed_file: dxpy.DXFile, fam_file: dxpy.DXFile, bim_file: dxpy.DXFile,
                             low_mac_list: dxpy.DXFile,
                             sparse_grm: dxpy.DXFile, sparse_grm_sample: dxpy.DXFile) -> None:
        # Now grab all genetic data that I have in the folder /project_resources/genetics/
        os.mkdir("genetics/")  # This is for legacy reasons to make sure all tests work...
        dxpy.download_dxfile(bed_file.get_id(), 'genetics/UKBB_470K_Autosomes_QCd.bed')
        dxpy.download_dxfile(bim_file.get_id(), 'genetics/UKBB_470K_Autosomes_QCd.bim')
        dxpy.download_dxfile(fam_file.get_id(), 'genetics/UKBB_470K_Autosomes_QCd.fam')
        dxpy.download_dxfile(low_mac_list.get_id(), 'genetics/UKBB_470K_Autosomes_QCd.low_MAC.snplist')
        # This is the sparse matrix
        dxpy.download_dxfile(sparse_grm.get_id(),
                             'genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx')
        dxpy.download_dxfile(sparse_grm_sample.get_id(),
                             'genetics/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt')

    @staticmethod
    def _process_regenie_snps(snp_qc_file: dxpy.DXFile) -> Optional[Path]:

        if snp_qc_file is None:
            return None
        else:
            dxpy.download_dxfile(snp_qc_file.get_id(), 'genetics/ukb_snp_qc.txt')
            with Path('genetics/ukb_snp_qc.txt').open('r') as snp_qc_reader,\
                    Path('genetics/rel_snps.txt').open('w') as rel_snps_writer:
                snp_qc_csv = csv.DictReader(snp_qc_reader, delimiter=" ")
                for snp in snp_qc_csv:
                    if snp['in_Relatedness'] == '1':
                        rel_snps_writer.write(f'{snp["rs_id"]}\n')
                snp_qc_reader.close()
                rel_snps_writer.close()

            return Path('genetics/rel_snps.txt')

    @staticmethod
    def _generate_filtered_genetic_data():

        # Generate a plink file to use that only has included individuals:
        cmd = "plink2 " \
              "--bfile /test/genetics/UKBB_470K_Autosomes_QCd --make-bed --keep-fam /test/SAMPLES_Include.txt " \
              "--out /test/genetics/UKBB_470K_Autosomes_QCd_WBA"
        run_cmd(cmd, True)

        # I have to do this to recover the sample information from plink
        cmd = "docker run -v /home/dnanexus/:/test/ egardner413/mrcepid-associationtesting plink2 " \
              "--bfile /test/genetics/UKBB_470K_Autosomes_QCd_WBA " \
              "--validate | grep samples"
        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()

        print(f'{"Plink individuals written":{65}}: {stdout.decode("utf-8").rstrip(" loaded from")}')
