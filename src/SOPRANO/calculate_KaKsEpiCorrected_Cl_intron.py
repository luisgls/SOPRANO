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

    site_names = ("NA", "NS")
    df_sites_extra = pd.read_csv(
        paths.epitope_nans, delimiter="\t", names=site_names
    )
    df_sites_intra = pd.read_csv(
        paths.intra_epitope_nans, delimiter="\t", names=site_names
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

    return merged_df, df_sites_extra, df_sites_intra


def _compute_mutation_counts(merged_df: pd.DataFrame):
    mutations_only = merged_df.drop(["EnsemblID"], axis=1)
    mutations_sums = mutations_only.sum(axis=0)
    return mutations_sums
