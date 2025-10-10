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


@pytest.mark.parametrize("tool, expected_output", [
    # ("regenie", "test.genes.REGENIE.stats.tsv"),
    # ("saige", "test.genes.SAIGE.stats.tsv"),
    ("bolt", "test.genes.BOLT.stats.tsv")
])
def test_burden_tools(tool, expected_output, temporary_path):
    """
    Parametrized test for burden tools (regenie, saige, bolt).
    This test will take a long time since it runs each tool entirely.
    :param tool: The tool to test (regenie, saige, bolt)
    :param expected_output: The expected output filename
    :param temporary_path: Temporary output path provided by pytest fixture
    """

    args = (
        f"--tool {tool} "
        f"--phenofile {test_data_dir}/phenotype.tsv "
        f"--association_tarballs {test_data_dir}/HC_PTV-MAF_001.tar.gz "
        f"--bgen_index {test_data_dir}/bgen_locs.tsv "
        f"--transcript_index {test_data_dir}/transcripts.tsv.gz "
        f"--base_covariates {test_data_dir}/base_covariates.covariates "
        f"--covarfile {test_data_dir}/other_covariates.covariates "
        f"--categorical_covariates batman "
        f"--array_bed_file {test_data_dir}/all_chroms.bed "
        f"--array_fam_file {test_data_dir}/all_chroms.fam "
        f"--array_bim_file {test_data_dir}/all_chroms.bim "
        f"--sparse_grm {test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx "
        f"--sparse_grm_sample {test_data_dir}/sparseGRM_470K_Autosomes_QCd.sparseGRM.mtx.sampleIDs.txt "
        f"--low_MAC_list {test_data_dir}/sim_chromosome.snplist "
    )

    if tool == "regenie":
        args += f"--regenie_run {test_data_dir}/step1_out2.tar.gz "

    load_class = LoadModule(output_prefix="test", input_args=args)
    load_class.start_module()

    expected_path = test_data_dir.parent / f"expected_outputs/{expected_output}"
    result_path = Path(expected_output)

    stats = pd.read_csv(expected_path, sep="\t")
    result = pd.read_csv(result_path, sep="\t")

    assert stats.shape == result.shape, f"Shape mismatch: {stats.shape} vs {result.shape}"
    assert stats.columns.tolist() == result.columns.tolist(), f"Columns mismatch: {stats.columns.tolist()} vs {result.columns.tolist()}"
    assert stats.equals(result), f"Values mismatch: {stats} vs {result}"
