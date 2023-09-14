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


def test__compute_mutation_counts():
    mock_df = pd.DataFrame(
        {
            "EnsemblID": ["ENSTxxx"] * 10,
            "intronrate": [123] * 10,
            "mut_1": [0] * 10,
            "mut_2": [1] * 10,
        }
    )

    expected_series = pd.Series(
        {"mut_1": 0, "mut_2": 10}, index=["mut_1", "mut_2"]
    )

    computed_series = dnds_intron._compute_mutation_counts(mock_df)

    assert computed_series.equals(expected_series)


def test__define_variables():
    mock_series = pd.Series({"m_1": 1, "m_2": 2})
    mock_df_1 = pd.DataFrame({"x_1": [3], "x_2": [4]})
    mock_df_2 = pd.DataFrame({"y_1": [5], "y_2": [6]})

    expected_series = pd.Series(
        {"m_1": 1, "m_2": 2, "x_1": 3, "x_2": 4, "y_1": 5, "y_2": 6}
    )

    computed_series = dnds_intron._define_variables(
        mock_series, mock_df_1, mock_df_2
    )

    assert computed_series.equals(expected_series)


def test__rescale_intron_by_synonymous():
    mock_vars = pd.Series(
        {
            "mutsintron": 1,
            "intra_synonymous_variant": 1,
            "extra_synonymous_variant": 1,
            "intra_site_2": 1,
            "extra_site_2": 1,
        }
    )

    expected_value = 1.0
    assert (
        dnds_intron._rescale_intron_by_synonymous(mock_vars) == expected_value
    )


def test__compute_kaks_intra_extra():
    mock_vars = pd.Series(
        {
            "intra_synonymous_variant": 1,
            "extra_synonymous_variant": 1,
            "intra_missense_variant": 2,
            "extra_missense_variant": 1,
            "intra_site_1": 1,
            "extra_site_1": 1,
            "intra_site_2": 2,
            "extra_site_2": 1,
        }
    )

    assert dnds_intron._compute_kaks_extra(mock_vars) == 1
    assert dnds_intron._compute_kaks_intra(mock_vars) == 4


def test__compute_kaks_intron():
    mock_vars = pd.Series(
        {
            "mutsintron": 1,
            "intra_synonymous_variant": 1,
            "extra_synonymous_variant": 1,
            "intra_missense_variant": 1,
            "extra_missense_variant": 1,
            "intra_site_1": 1,
            "extra_site_1": 1,
            "intra_site_2": 1,
            "extra_site_2": 1,
        }
    )

    assert dnds_intron._compute_kaks_intron(mock_vars) == 1.0


def test__compute_conf_interval():
    # TODO: Implement
    pass


merged, extra, intra = dnds_intron._preprocess_dfs(paths)

muts = dnds_intron._compute_mutation_counts(merged)

vars = dnds_intron._define_variables(muts, extra, intra)

cl_epi = dnds_intron._compute_conf_interval(vars, "extra", "katz")
cl_nonepi = dnds_intron._compute_conf_interval(vars, "intra", "katz")
cl_intron = dnds_intron._compute_conf_interval(vars, "intron", "katz")

pval_nonepi = dnds_intron._compute_pvalue(vars, cl_nonepi, "intra")
pval_intron = dnds_intron._compute_pvalue(vars, cl_nonepi, "intron")

print("epitope", cl_epi)
print("non-epitope", cl_nonepi)
print("intron", cl_intron)
print("pval non-epitope", pval_nonepi)
print("pval intron", pval_intron)
