import pathlib

import pandas as pd

from SOPRANO import calculate_KaKsEpiCorrected_Cl_intron as dnds_intron
from SOPRANO.objects import AnalysisPaths

# TODO: Fixture....
data_dir = pathlib.Path("dnds_data")
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


def test__preprocess_dfs():
    merged, sites_extra, sites_intra = dnds_intron._preprocess_dfs(paths)

    test_merged_path = data_dir.joinpath("merged.tsv")
    test_extra_path = data_dir.joinpath("sites_extra.tsv")
    test_intra_path = data_dir.joinpath("sites_intra.tsv")

    assert merged.head().equals(pd.read_csv(test_merged_path, delimiter="\t"))
    assert sites_extra.head().equals(
        pd.read_csv(test_extra_path, delimiter="\t")
    )
    assert sites_intra.head().equals(
        pd.read_csv(test_intra_path, delimiter="\t")
    )
