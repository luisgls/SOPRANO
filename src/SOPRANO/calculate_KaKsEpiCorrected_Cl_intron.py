import numpy as np
import pandas as pd

from SOPRANO.objects import AnalysisPaths

_DISCARD_INTRON = True

# TODO: Discuss - this looks like a potential bug
#       in the R implementation...
#       df3<-merge(df2,df.intron,by.x="EnsembleID",by.y="EnsemblID", all.x = T)
#       Effectively, this removes the intron values from the calculation.
#       Probably want all.y=T where y is df.intron
#       The above flag mimics this behaviour.


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

    if _DISCARD_INTRON:
        merged_df["mutsintron"] *= 0
        merged_df["intronlength"] *= 0

    return merged_df, sites_extra_df, sites_intra_df


def _compute_mutation_counts(merged_df: pd.DataFrame) -> pd.Series:
    mutations_only = merged_df.drop(["EnsemblID", "intronrate"], axis=1)
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


def _rescale_intron_by_synonymous(variables: pd.Series):
    # Sum over intron and synonymous muts
    n_intron = variables["mutsintron"]
    n_intra_syn = variables["intra_synonymous_variant"]
    n_extra_syn = variables["extra_synonymous_variant"]

    # Second element of intra and extra sites
    n_intra_sites_2 = variables["intra_site_2"]
    n_extra_sites_2 = variables["extra_site_2"]

    ri = n_intra_syn / n_intra_sites_2
    re = n_extra_syn / n_extra_sites_2

    return 2 * n_intron / (ri + re)


def _compute_kaks(variables: pd.Series, prefix: str):
    if prefix not in ("intra", "extra"):
        raise ValueError(prefix)

    # Sum vals over missesne and synonymous muts
    n_mis = variables[f"{prefix}_missense_variant"]
    n_syn = variables[f"{prefix}_synonymous_variant"]

    # First and second elements of sites vector
    sites_1 = variables[f"{prefix}_site_1"]
    sites_2 = variables[f"{prefix}_site_2"]

    rn = n_mis / n_syn
    rs = sites_2 / sites_1

    return rn * rs


def _compute_kaks_extra(variables: pd.Series):
    return _compute_kaks(variables, prefix="extra")


def _compute_kaks_intra(variables: pd.Series):
    return _compute_kaks(variables, prefix="intra")


def _compute_kaks_intron(variables: pd.Series):
    # Sum vals over intron, missesne and synonymous intra muts
    n_int = variables["mutsintron"]
    n_mis = variables["intra_missense_variant"]
    n_syn = variables["intra_synonymous_variant"]

    # First and second elements of sites vector
    sites_1 = variables["intra_site_1"]
    sites_2 = variables["intra_site_2"]

    # Get intron muts rescaled by synonymous muts
    n_intron_rescaled = _rescale_intron_by_synonymous(variables)

    rn = n_mis / (n_syn + n_int)
    rs = (sites_2 + n_intron_rescaled) / sites_1

    return rn * rs


def _katz_confidence_interval(n_mis, n_syn, sites_1, sites_2, prefix):
    # TODO: Find out why +1 in sites:
    # Is this to avoid divergence in case of 0?
    # If so, maybe consider max(1, sites_x)
    p1 = n_mis / (sites_1 + 1)
    p2 = n_syn / (sites_2 + 1)

    global_dnds = p1 / p2

    f1 = (1 - p1) / (sites_1 * p1)
    f2 = (1 - p2) / (sites_2 * p2)

    sqrt = np.sqrt(f1 + f2)

    low = global_dnds * np.exp(-1.96 * sqrt)
    high = global_dnds * np.exp(1.96 * sqrt)

    return pd.Series({f"{prefix}_Cl_low": low, f"{prefix}_Cl_high": high})


def _compute_conf_interval(variables: pd.Series, prefix: str, method: str):
    if method != "katz":
        raise ValueError(f"Unimplemented method: {method}")

    if prefix == "intra":
        n_mis = variables["intra_missense_variant"]
        n_syn = variables["intra_synonymous_variant"]
        sites_1 = variables["intra_site_1"]
        sites_2 = variables["intra_site_2"]
    elif prefix == "extra":
        n_mis = variables["extra_missense_variant"]
        n_syn = variables["extra_synonymous_variant"]
        sites_1 = variables["extra_site_1"]
        sites_2 = variables["extra_site_2"]
    elif prefix == "intron":
        n_mis = variables["intra_missense_variant"]
        n_syn = variables["intra_synonymous_variant"] + variables["mutsintron"]
        sites_1 = variables["intra_site_1"]
        sites_2 = variables["intra_site_2"] + _rescale_intron_by_synonymous(
            variables
        )
    else:
        raise ValueError(prefix)

    return _katz_confidence_interval(n_mis, n_syn, sites_1, sites_2, prefix)


def _compute_pvalue(
    variables: pd.Series, conf_intervals: pd.Series, prefix: str
):
    """

    :param variables:
    :param conf_intervals: Confidence intervals computed from
            non-target region with or without intronic information
    """
    high = conf_intervals["intra_Cl_high"]
    low = conf_intervals["intra_Cl_low"]

    extra_kak = _compute_kaks_extra(variables)

    if prefix == "intron":
        intra_kak = _compute_kaks_intron(variables)
    elif prefix == "intra":
        intra_kak = _compute_kaks_intra(variables)
    else:
        raise ValueError(f"pvalue not configured for prefix: {prefix}")

    delta_kak = extra_kak - intra_kak
    standard_error = (high - low) / 3.92
    z = delta_kak / standard_error

    pval_1 = np.exp(-0.717 * z - 0.416 * z**2)
    pval_2 = np.exp(0.717 * z + 0.416 * z**2)

    if 0 < pval_2 <= 1:
        pval = pval_2
    else:
        pval = pval_1

    pval = max([pval, 1e-4])

    return pval
