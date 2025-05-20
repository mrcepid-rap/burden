"""
This runs the burden association tests using the LoadModule class. Note that if you are running REGENIE step1 for
these tests, they will take a while. Once you have run step1 once, you can save the output files to save time
on running step2. Step1 files are also provided in case you want to copy them into the test directory and avoid
running step1.

Note that BOLT will always take a long time to run.

For access to the test data, please contact one of the authors of the burden package.
"""

import shutil
from pathlib import Path

import pandas as pd
import pytest

from burden.loader import LoadModule

# Set this flag to True if you want to keep (copy) the temporary output files
KEEP_TEMP = True

# test data directory
test_data_dir = Path(__file__).parent / 'test_data'


@pytest.fixture
def temporary_path(tmp_path, monkeypatch):
    """
    Prepare a temporary working directory that contains a copy of the test_data
    directory, then change the working directory to it.

    If KEEP_TEMP is True, after the test the entire temporary directory will be copied
    to a folder 'temp_test_outputs' in the project root.
    """
    # Determine where the original test_data directory is located.
    # (Assumes it is at <project_root>/test_data)
    test_data_source = Path(__file__).parent / "test_data"

    # Create the destination folder inside the tmp_path.
    destination = tmp_path / "test_data"
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Copy the entire test_data directory into the temporary directory.
    shutil.copytree(test_data_source, destination)

    # Change the current working directory to the temporary directory.
    monkeypatch.chdir(tmp_path)

    # Yield the temporary directory to the test.
    yield tmp_path

    # After the test, if KEEP_TEMP is True, copy the temporary directory to a persistent location.
    if KEEP_TEMP:
        persistent_dir = Path(__file__).parent / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(exist_ok=True)
        shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)
        print(f"Temporary output files have been copied to: {persistent_dir}")


def test_regenie(temporary_path) -> None:
    """
    Test the LoadModule class with a set of input arguments.

    Note that this test will take a very long time to run, as it is running the runners (bolt, regenie etc.) entirely.
    You may wish to save the output of regenie by setting KEEP_TEMP=True, so that you can save the output of regenie step1
    to your /test_data/ directory (make sure to add it to .gitignore). This will greatly speed up the time to run regenie for
    testing purposes.

    :param temporary_path: where the temporary output files are going to be stored
    :return: None
    """

    # if regenie step1 files exist, copy them into our current pytest temp directory
    # see if they exist
    if (Path(test_data_dir / 'fit_out_pred.list').exists() and Path(
            test_data_dir / 'fit_out_pred.list').stat().st_size > 0) and \
            (Path(test_data_dir / 'fit_out_1.loco').exists() and Path(
                test_data_dir / 'fit_out_1.loco').stat().st_size > 0):
        # copy the files from the test_data directory to the current pytest directory
        shutil.copy(test_data_dir / 'fit_out_pred.list', 'fit_out_pred.list')
        shutil.copy(test_data_dir / 'fit_out_1.loco', 'fit_out_1.loco')

    # run the burden class
    load_class = LoadModule(
        output_prefix="test",
        input_args=(
            "--tool "
            "regenie "
            "--phenofile "
            f"{test_data_dir}/phenotype.tsv "
            "--association_tarballs "
            f"{test_data_dir}/HC_PTV-MAF_001.tar.gz "
            "--bgen_index "
            f"{test_data_dir}/bgen_locs.tsv "
            # "--regenie_smaller_snps "
            # "/Users/alish.palmos/PycharmProjects/burden/test/test_data/ukb_snp_qc.txt "
            "--transcript_index "
            f"{test_data_dir}/transcripts.tsv.gz "
            "--base_covariates "
            f"{test_data_dir}/base_covariates.covariates "
            "--covarfile "
            f"{test_data_dir}/other_covariates.covariates "
            "--categorical_covariates "
            "batman "
            "--array_bed_file "
            f"{test_data_dir}/all_chroms.bed "
            "--array_fam_file "
            f"{test_data_dir}/all_chroms.fam "
            "--array_bim_file "
            f"{test_data_dir}/all_chroms.bim "
            "--sparse_grm "
            f"{test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx "
            "--sparse_grm_sample "
            f"{test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt "
            "--low_MAC_list "
            f"{test_data_dir}/sim_chromosome.snplist "
        )
    )

    load_class.start_module()

    # check the output files
    # read in expected
    stats = pd.read_csv(test_data_dir.parent / "expected_outputs/test.genes.REGENIE.stats.tsv", sep="\t")
    # read in the output
    result = pd.read_csv("test.genes.REGENIE.stats.tsv", sep="\t")
    # check the output
    assert stats.shape == result.shape, f"Shape mismatch: {stats.shape} vs {result.shape}"
    # check the columns
    assert stats.columns.tolist() == result.columns.tolist(), f"Columns mismatch: {stats.columns.tolist()} vs {result.columns.tolist()}"
    # check the values
    assert stats.equals(result), f"Values mismatch: {stats} vs {result}"


def test_saige(temporary_path) -> None:
    """
    Test the LoadModule class with a set of input arguments.

    Note that this test will take a very long time to run, as it is running the runners (bolt, regenie etc.) entirely.
    You may wish to save the output of regenie by setting KEEP_TEMP=True, so that you can save the output of regenie step1
    to your /test_data/ directory (make sure to add it to .gitignore). This will greatly speed up the time to run regenie for
    testing purposes.

    :param temporary_path: where the temporary output files are going to be stored
    :return: None
    """

    # run the burden class
    load_class = LoadModule(
        output_prefix="test",
        input_args=(
            "--tool "
            "saige "
            "--phenofile "
            f"{test_data_dir}/phenotype.tsv "
            "--association_tarballs "
            f"{test_data_dir}/HC_PTV-MAF_001.tar.gz "
            "--bgen_index "
            f"{test_data_dir}/bgen_locs.tsv "
            # "--regenie_smaller_snps "
            # "/Users/alish.palmos/PycharmProjects/burden/test/test_data/ukb_snp_qc.txt "
            "--transcript_index "
            f"{test_data_dir}/transcripts.tsv.gz "
            "--base_covariates "
            f"{test_data_dir}/base_covariates.covariates "
            "--covarfile "
            f"{test_data_dir}/other_covariates.covariates "
            "--categorical_covariates "
            "batman "
            "--array_bed_file "
            f"{test_data_dir}/all_chroms.bed "
            "--array_fam_file "
            f"{test_data_dir}/all_chroms.fam "
            "--array_bim_file "
            f"{test_data_dir}/all_chroms.bim "
            "--sparse_grm "
            f"{test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx "
            "--sparse_grm_sample "
            f"{test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt "
            "--low_MAC_list "
            f"{test_data_dir}/sim_chromosome.snplist "
        )
    )

    load_class.start_module()

    # check the output files
    # read in expected
    stats = pd.read_csv(test_data_dir.parent / "expected_outputs/test.genes.SAIGE.stats.tsv", sep="\t")
    # read in the output
    result = pd.read_csv("test.genes.SAIGE.stats.tsv", sep="\t")
    # check the output
    assert stats.shape == result.shape, f"Shape mismatch: {stats.shape} vs {result.shape}"
    # check the columns
    assert stats.columns.tolist() == result.columns.tolist(), f"Columns mismatch: {stats.columns.tolist()} vs {result.columns.tolist()}"
    # check the values
    assert stats.equals(result), f"Values mismatch: {stats} vs {result}"


def test_bolt(temporary_path) -> None:
    """
    Test the LoadModule class with a set of input arguments.

    Note that this test will take a very long time to run, as it is running the runners (bolt, regenie etc.) entirely.
    You may wish to save the output of regenie by setting KEEP_TEMP=True, so that you can save the output of regenie step1
    to your /test_data/ directory (make sure to add it to .gitignore). This will greatly speed up the time to run regenie for
    testing purposes.

    :param temporary_path: where the temporary output files are going to be stored
    :return: None
    """

    # run the burden class
    load_class = LoadModule(
        output_prefix="test",
        input_args=(
            "--tool "
            "bolt "
            "--phenofile "
            f"{test_data_dir}/phenotype.tsv "
            "--association_tarballs "
            f"{test_data_dir}/HC_PTV-MAF_001.tar.gz "
            "--bgen_index "
            f"{test_data_dir}/bgen_locs.tsv "
            # "--regenie_smaller_snps "
            # "/Users/alish.palmos/PycharmProjects/burden/test/test_data/ukb_snp_qc.txt "
            "--transcript_index "
            f"{test_data_dir}/transcripts.tsv.gz "
            "--base_covariates "
            f"{test_data_dir}/base_covariates.covariates "
            "--covarfile "
            f"{test_data_dir}/other_covariates.covariates "
            "--categorical_covariates "
            "batman "
            "--array_bed_file "
            f"{test_data_dir}/all_chroms.bed "
            "--array_fam_file "
            f"{test_data_dir}/all_chroms.fam "
            "--array_bim_file "
            f"{test_data_dir}/all_chroms.bim "
            "--sparse_grm "
            f"{test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx "
            "--sparse_grm_sample "
            f"{test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt "
            "--low_MAC_list "
            f"{test_data_dir}/sim_chromosome.snplist "
        )
    )

    load_class.start_module()

    # check the output files
    # read in expected
    stats = pd.read_csv(test_data_dir.parent / "expected_outputs/test.genes.BOLT.stats.tsv", sep="\t")
    # read in the output
    result = pd.read_csv("test.genes.BOLT.stats.tsv", sep="\t")
    # check the output
    assert stats.shape == result.shape, f"Shape mismatch: {stats.shape} vs {result.shape}"
    # check the columns
    assert stats.columns.tolist() == result.columns.tolist(), f"Columns mismatch: {stats.columns.tolist()} vs {result.columns.tolist()}"
    # check the values
    assert stats.equals(result), f"Values mismatch: {stats} vs {result}"
