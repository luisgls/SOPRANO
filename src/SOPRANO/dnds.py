import numpy as np
import pandas as pd

from SOPRANO.objects import (
    AnalysisPaths,
    AuxiliaryFiles,
    AuxiliaryPaths,
    Parameters,
)
from SOPRANO.pipeline_utils import (
    MissingDataError,
    _PipelineComponent,
    is_empty,
)
from SOPRANO.sh_utils import subprocess_pipes


def _intersect_introns(paths: AnalysisPaths, aux_files: AuxiliaryPaths):
    """
    Implements:

    intersectBed -a $BASEDIR/data/transcript_intron_length.bed -b
    $TMP/$NAME.intronic.bed -wo |
        mergeBed -i stdin -c 4,5,6,10,11 -o mode,mode,mode,collapse,count |
            awk '{print $4"\t"$8/($6+1)"\t"$8"\t"$6}' >
                $TMP/$NAME.intronic.rate

    :param paths:
    :return:
    """

    subprocess_pipes.pipe(
        [
            "intersectBed",
            "-a",
            aux_files.intron_length.as_posix(),
            "-b",
            paths.variants_intronic.as_posix(),
            "-wo",
        ],
        [
            "mergeBed",
            "-i",
            "stdin",
            "-c",
            "4,5,6,10,11",
            "-o",
            "mode,mode,mode,collapse,count",
        ],
        ["awk", r'{print $4"\t"$8/($6+1)"\t"$8"\t"$6}'],
        output_path=paths.intron_rate,
    )


class ComputeIntronRate(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        if not params.variants_intronic.exists():
            raise MissingDataError(params.variants_intronic)

    @staticmethod
    def apply(params: Parameters):
        _intersect_introns(params, AuxiliaryFiles)


class ComputeStatistics(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.data_epitopes,
            params.epitope_nans,
            params.intra_epitope_nans,
            params.intron_rate,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ComputeStatistics.check_ready(params)
        _compute_coverage(params)


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

    if is_empty(paths.intron_rate):
        merged_df = spread_mutations_df
    else:
        intron_df = pd.read_csv(
            paths.intron_rate,
            delimiter="\t",
            names=("EnsemblID", "intronrate", "mutsintron", "intronlength"),
        )

        # Merge mutations with intron df
        merged_df = spread_mutations_df.merge(
            intron_df, how="outer", on="EnsemblID"
        ).fillna(0)

        # TODO: Double check this is equivalent to
        #   merge(df2,df.intron,by.x="EnsembleID",by.y="EnsemblID", all.x = T)

    return merged_df, sites_extra_df, sites_intra_df


def _compute_mutation_counts(merged_df: pd.DataFrame) -> pd.Series:
    try:
        mutations_only = merged_df.drop(["EnsemblID", "intronrate"], axis=1)
    except KeyError:
        mutations_only = merged_df.drop(["EnsemblID"], axis=1)
    mutations_sums = mutations_only.sum(axis=0)  # type: pd.Series

    n_total_epitope = (
        mutations_sums["extra_missense_variant"]
        + mutations_sums["extra_synonymous_variant"]
    )
    n_total_nonepitope = (
        mutations_sums["intra_missense_variant"]
        + mutations_sums["intra_synonymous_variant"]
    )
    combined_sums = pd.Series(
        {
            "mut_total_epitope": n_total_epitope,
            "mut_total_non_epitope": n_total_nonepitope,
        }
    )
    if "mutsintron" in mutations_sums:
        n_total_intron = n_total_nonepitope + mutations_sums["mutsintron"]
        combined_sums = pd.concat(
            [combined_sums, pd.Series({"mut_total_intron": n_total_intron})]
        )

    return pd.concat([mutations_sums, combined_sums])


def _define_variables(
    mutation_sums: pd.Series,
    extra_sites: pd.DataFrame,
    intra_sites: pd.DataFrame,
) -> pd.Series:
    variables = mutation_sums.copy()  # type: pd.Series
    sq_extra = extra_sites.squeeze()  # type: pd.Series
    sq_intra = intra_sites.squeeze()  # type: pd.Series
    variables = pd.concat([variables, sq_extra])
    variables = pd.concat([variables, sq_intra])
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
    # Usually sites number is >> 1, but it is possible to have 0 value.
    # Therefore, take max(sites, 1) to avoid divergence in ratios
    sites_1 = max([sites_1, 1])
    sites_2 = max([sites_2, 1])

    p1 = n_mis / sites_1
    p2 = n_syn / sites_2

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


def _compute_all_conf_intervals(variables: pd.Series, method: str):
    conf_intervals_extra_intra = pd.concat(
        [
            _compute_conf_interval(variables, "extra", method),
            _compute_conf_interval(variables, "intra", method),
        ]
    )

    if "mutsintron" in variables:
        return pd.concat(
            [
                conf_intervals_extra_intra,
                _compute_conf_interval(variables, "intron", method),
            ]
        )
    else:
        return conf_intervals_extra_intra


def _compute_pvalue(
    variables: pd.Series, conf_intervals: pd.Series, prefix: str
):
    """

    :param variables:
    :param conf_intervals: Confidence intervals computed from
            non-target region with or without intronic information
    """

    extra_kak = _compute_kaks_extra(variables)

    if prefix == "intron":
        intra_kak = _compute_kaks_intron(variables)
        high = conf_intervals["intron_Cl_high"]
        low = conf_intervals["intron_Cl_low"]
    elif prefix == "intra":
        intra_kak = _compute_kaks_intra(variables)
        high = conf_intervals["intra_Cl_high"]
        low = conf_intervals["intra_Cl_low"]
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


def _compute_coverage(paths: AnalysisPaths):
    merged_df, extra_df, intra_df = _preprocess_dfs(paths)

    mut_counts = _compute_mutation_counts(merged_df)

    vars = _define_variables(mut_counts, extra_df, intra_df)

    kaks_extra = _compute_kaks_extra(vars)
    kaks_intra = _compute_kaks_intra(vars)

    conf_intervals = _compute_all_conf_intervals(vars, "katz")
    pval_intra = _compute_pvalue(vars, conf_intervals, "intra")

    results_df = pd.DataFrame(
        {
            "Coverage": ["Exonic_Only"],
            "ON_dNdS": [kaks_extra],
            "ON_Low_CI": [conf_intervals["extra_Cl_low"]],
            "ON_High_CI": [conf_intervals["extra_Cl_high"]],
            "ON_Mutations": [mut_counts["mut_total_epitope"]],
            "OFF_dNdS": [kaks_intra],
            "OFF_Low_CI": [conf_intervals["intra_Cl_low"]],
            "OFF_High_CI": [conf_intervals["intra_Cl_high"]],
            "OFF_Mutations": [mut_counts["mut_total_non_epitope"]],
            "Pvalue": [pval_intra],
            "ON_na": [vars["extra_missense_variant"]],
            "ON_NA": [vars["extra_site_1"]],
            "ON_ns": [vars["extra_synonymous_variant"]],
            "ON_NS": [vars["extra_site_2"]],
            "OFF_na": [vars["intra_missense_variant"]],
            "OFF_NA": [vars["intra_site_1"]],
            "OFF_ns": [vars["intra_synonymous_variant"]],
            "OFF_NS": [vars["intra_site_2"]],
        }
    )

    if "mutsintron" in mut_counts:
        kaks_intron = _compute_kaks_intron(vars)
        pval_intron = _compute_pvalue(vars, conf_intervals, "intron")

        rs_intron = _rescale_intron_by_synonymous(vars)

        intron_df = pd.DataFrame(
            {
                "Coverage": ["Exonic_Intronic"],
                "ON_dNdS": [kaks_extra],
                "ON_Low_CI": [conf_intervals["extra_Cl_low"]],
                "ON_High_CI": [conf_intervals["extra_Cl_high"]],
                "ON_Mutations": [mut_counts["mut_total_epitope"]],
                "OFF_dNdS": [kaks_intron],
                "OFF_Low_CI": [conf_intervals["intron_Cl_low"]],
                "OFF_High_CI": [conf_intervals["intron_Cl_high"]],
                "OFF_Mutations": [mut_counts["mut_total_intron"]],
                "Pvalue": [pval_intron],
                "ON_na": [vars["extra_missense_variant"]],
                "ON_NA": [vars["extra_site_1"]],
                "ON_ns": [vars["extra_synonymous_variant"]],
                "ON_NS": [vars["extra_site_2"]],
                "OFF_na": [vars["intra_missense_variant"]],
                "OFF_NA": [vars["intra_site_1"]],
                "OFF_ns": [
                    vars["intra_synonymous_variant"] + vars["mutsintron"]
                ],
                "OFF_NS": [vars["intra_site_2"] + rs_intron],
            }
        )
        results_df = pd.concat([results_df, intron_df], axis=0)

        # NOTE: ON/OFF_NS appear sensitive to the reference genome release!

    print(f"Exporting results to {paths.results_path}:")
    print(results_df)
    results_df.to_csv(paths.results_path, sep="\t")
