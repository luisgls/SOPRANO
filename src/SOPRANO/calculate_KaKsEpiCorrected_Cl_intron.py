import pandas as pd

from SOPRANO.objects import AnalysisPaths


def _preprocess_dfs(paths: AnalysisPaths):
    data_path = paths.data_epitopes.as_posix()

    # Load file as data frame and assign headers
    data_epitopes_df = pd.read_csv(
        data_path, delimiter="\t", names=("EnsemblID", "Total", "Class")
    )

    # Spread possible values in the mutation classes to have their own
    # columns. Since these numbers are mutually exclusive,
    # i.e., can't have missense and extra synonymous with non-zero values,
    # we fill the missing data with zeros.
    spread_mutations_df = (
        data_epitopes_df.pivot(
            index="EnsemblID", columns="Class", values="Total"
        )
        .fillna(0)
        .reset_index()
    )

    site_names = ("site_1", "site_2")
    extra_site_names = ["extra_" + s for s in site_names]
    intra_site_names = ["intra_" + s for s in site_names]

    sites_extra_df = pd.read_csv(
        paths.epitope_nans, delimiter="\t", names=extra_site_names
    )
    sites_intra_df = pd.read_csv(
        paths.intra_epitope_nans, delimiter="\t", names=intra_site_names
    )

    intron_df = pd.read_csv(
        paths.intron_rate,
        delimiter="\t",
        names=("EnsemblID", "intronrate", "mutsintron", "intronlength"),
    )

    # Merge mutations with intron df
    merged_df = spread_mutations_df.merge(
        intron_df, how="outer", on="EnsemblID"
    ).fillna(0)

    return merged_df, sites_extra_df, sites_intra_df


def _compute_mutation_counts(merged_df: pd.DataFrame) -> pd.Series:
    mutations_only = merged_df.drop(["EnsemblID"], axis=1)
    mutations_sums = mutations_only.sum(axis=0)
    return mutations_sums


def _define_variables(
    mutation_sums: pd.Series,
    extra_sites: pd.DataFrame,
    intra_sites: pd.DataFrame,
):
    variables = mutation_sums.copy()

    variables = pd.concat([variables, extra_sites.squeeze()])
    variables = pd.concat([variables, intra_sites.squeeze()])
    return variables
