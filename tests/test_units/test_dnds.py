import pandas as pd

import SOPRANO.core.dnds


def test__preprocess_dfs(tcga_05_4396_ssb192_cfg):
    (
        paths,
        merged_path,
        sites_extra_path,
        sites_intra_path,
    ) = tcga_05_4396_ssb192_cfg

    merged, sites_extra, sites_intra = SOPRANO.core.dnds._preprocess_dfs(paths)

    assert merged.head().equals(pd.read_csv(merged_path, delimiter="\t"))
    assert sites_extra.head().equals(
        pd.read_csv(sites_extra_path, delimiter="\t")
    )
    assert sites_intra.head().equals(
        pd.read_csv(sites_intra_path, delimiter="\t")
    )


def test__compute_mutation_counts():
    mock_df = pd.DataFrame(
        {
            "EnsemblID": ["ENSTxxx"] * 10,
            "intronrate": [123] * 10,
            "extra_missense_variant": [0] * 10,
            "extra_synonymous_variant": [1] * 10,
            "intra_missense_variant": [0] * 10,
            "intra_synonymous_variant": [1] * 10,
        }
    )

    expected_series = pd.Series(
        {
            "extra_missense_variant": 0,
            "extra_synonymous_variant": 10,
            "intra_missense_variant": 0,
            "intra_synonymous_variant": 10,
            "mut_total_epitope": 10,
            "mut_total_non_epitope": 10,
        },
        index=[
            "extra_missense_variant",
            "extra_synonymous_variant",
            "intra_missense_variant",
            "intra_synonymous_variant",
            "mut_total_epitope",
            "mut_total_non_epitope",
        ],
    )

    computed_series = SOPRANO.core.dnds._compute_mutation_counts(mock_df)

    assert computed_series.equals(expected_series)


def test__define_variables():
    mock_series = pd.Series({"m_1": 1, "m_2": 2})
    mock_df_1 = pd.DataFrame({"x_1": [3], "x_2": [4]})
    mock_df_2 = pd.DataFrame({"y_1": [5], "y_2": [6]})

    expected_series = pd.Series(
        {"m_1": 1, "m_2": 2, "x_1": 3, "x_2": 4, "y_1": 5, "y_2": 6}
    )

    computed_series = SOPRANO.core.dnds._define_variables(
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
        SOPRANO.core.dnds._rescale_intron_by_synonymous(mock_vars)
        == expected_value
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

    assert SOPRANO.core.dnds._compute_kaks_extra(mock_vars) == 1
    assert SOPRANO.core.dnds._compute_kaks_intra(mock_vars) == 4


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

    assert SOPRANO.core.dnds._compute_kaks_intron(mock_vars) == 1.0


def test__compute_conf_interval():
    # TODO: Implement
    pass
