import pathlib

import pytest

from SOPRANO.core.objects import AnalysisPaths


@pytest.fixture
def tcga_05_4396_ssb192_cfg():
    data_dir = pathlib.Path(__file__).parent.joinpath("dnds_data")
    data_epitopes_path = data_dir.joinpath("TCGA-05-4396.data.epitopes")
    epi_nans_path = data_dir.joinpath("TCGA-05-4396.epitopes.nans")
    intra_nans_path = data_dir.joinpath("TCGA-05-4396.intra_epitopes.nans")
    intron_path = data_dir.joinpath("TCGA-05-4396.intron.rate")

    merged_path = data_dir.joinpath("merged.tsv")
    sites_extra_path = data_dir.joinpath("sites_extra.tsv")
    sites_intra_path = data_dir.joinpath("sites_intra.tsv")

    paths = AnalysisPaths(
        "tcga_test", data_epitopes_path, data_epitopes_path, pathlib.Path(".")
    )  # Just hacking something together here...

    paths.data_epitopes = data_epitopes_path
    paths.epitope_nans = epi_nans_path
    paths.intra_epitope_nans = intra_nans_path
    paths.intron_rate = intron_path

    return paths, merged_path, sites_extra_path, sites_intra_path
