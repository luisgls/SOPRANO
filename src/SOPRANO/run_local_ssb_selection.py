import argparse
import pathlib
import subprocess
from datetime import datetime

from SOPRANO import (
    analysis,
    intersect,
    objects,
    obtain_fasta_regions,
    prepare_coordinates,
)


def check_path(cli_path: pathlib.Path | None, optional=False):
    if cli_path is None:
        if not optional:
            raise Exception("Input path is not optional and path is None!")
    elif not cli_path.exists():
        raise Exception(f"CLI input path does not exist: {cli_path}")


def check_genome(_namespace: argparse.Namespace) -> str:
    ref = _namespace.genome_ref

    available_refs = ("grch37", "grch38")

    if ref.lower() not in available_refs:
        raise ValueError(
            f"Reference {ref} not supported. "
            f"Permitted choices: {available_refs}"
        )

    return ref.lower()


def parse_args():
    parser = argparse.ArgumentParser(description="SOPRANO input arguments")

    parser.add_argument(
        "--input",
        "-i",
        dest="input_path",
        type=pathlib.Path,
        help="Prove the path to the input VEP annotated file.",
        required=True,
    )

    parser.add_argument(
        "--bed_file",
        "-b",
        dest="bed_path",
        type=pathlib.Path,
        help="Provide the path to the bed file with protein coordinates named "
        "by Transcript (ENSTXXXXXX 123 135)",
        required=True,
    )

    parser.add_argument(
        "--output",
        "-o",
        dest="cache_dir",
        type=pathlib.Path,
        help="Provide the path to the output directory in which dN/dS results "
        "will be cached.",
        required=True,
    )

    parser.add_argument(
        "--name",
        "-n",
        dest="analysis_name",
        type=str,
        help="Provide an identifying name for your results.",
        required=True,
    )
    analysis_params_group = parser.add_argument_group()

    analysis_params_group.add_argument(
        "--target_regions",
        "-t",
        dest="target_regions",
        type=pathlib.Path,
        help="Provide a bed file with regions to randomize.",
    )

    analysis_params_group.add_argument(
        "--use_ssb192",
        dest="use_ssb192",
        action="store_true",
        help="If flag is used ssb192 will be used, otherwise ssb7.",
    )

    analysis_params_group.add_argument(
        "--use_random",
        dest="use_random",
        action="store_true",
        help="If flag is used, calculate a dNdS calue for a random region "
        "similar to the target.",
    )

    analysis_params_group.add_argument(
        "--exclude_drivers",
        dest="exclude_drivers",
        action="store_true",
        help="If flag is used, driver geners will be excluded from the "
        "calculation.",
    )

    analysis_params_group.add_argument(
        "--seed",
        dest="seed",
        default=-1,
        type=int,
        help="Provide seed value for shuffle process in randomization. By "
        "default, seed value is < 0, for which, no seed value will be "
        "applied.",
    )

    transcript_args = parser.add_argument_group()

    transcript_args.add_argument(
        "--transcript",
        dest="transcript",
        help="Provide path to transcript file",
        default=objects.EnsemblTranscripts.transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--protein_transcript",
        dest="protein_transcript",
        help="Provide path to protein transcript file",
        default=objects.EnsemblTranscripts.protein_transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--fasta",
        dest="transcript_ids",
        help="Provide path to the ensembl transcript IDs fasta file",
        default=objects.EnsemblTranscripts.transcript_fasta,
        type=pathlib.Path,
    )

    genome_args = parser.add_argument_group()

    genome_args.add_argument(
        "--reference",
        "-r",
        dest="genome_ref",
        default="GRCh37",
        type=str,
        help="Reference genome file definition. By default, uses GRCh37. Pass "
        "instead GRCh38 if preferred. In order to download the reference "
        "fasta file, you can execute the command:\n"
        "GET_GENOMES -r GRCh37\n"
        "or\n"
        "GET_GENOMES -r GRCh38",
    )

    args = parser.parse_args()

    check_path(args.input_path)
    check_path(args.bed_path)
    check_path(args.cache_dir)
    check_path(args.target_regions, optional=True)
    check_path(args.transcript)
    check_path(args.protein_transcript)
    check_path(args.transcript_ids)
    check_genome(args)

    return args


def title_output():
    print(
        """
███████  ██████  ██████  ██████   █████  ███    ██  ██████  
██      ██    ██ ██   ██ ██   ██ ██   ██ ████   ██ ██    ██ 
███████ ██    ██ ██████  ██████  ███████ ██ ██  ██ ██    ██ 
     ██ ██    ██ ██      ██   ██ ██   ██ ██  ██ ██ ██    ██ 
███████  ██████  ██      ██   ██ ██   ██ ██   ████  ██████  
"""
    )


def time_output():
    now = datetime.now()
    return now.strftime("%d/%m/%Y %H:%M:%S")


def task_output(msg):
    print(f"[{time_output()}] {msg}")


def line_output(n=60):
    print("-" * n)


def startup_output(**kwargs):
    title_output()
    line_output()
    print("Selection On PRotein ANnotated regiOns")
    line_output()

    if kwargs:
        # Parameters used in pipeline
        for k, v in kwargs.items():
            print("-> {0:.<30}".format(k) + f"{v}")

        line_output()


def main(_namespace=None):
    cli_args = parse_args() if _namespace is None else _namespace
    startup_output(**cli_args.__dict__)
    params = objects.Parameters.from_namespace(cli_args)

    task_output("Filtering transcripts")
    prepare_coordinates.FilterTranscripts.apply(params)

    if params.use_target_regions:
        task_output(
            f"Randomizing transcripts subject to bed file: "
            f"{params.target_regions_path.as_posix()}"
        )
        randomization_method = prepare_coordinates.RandomizeWithRegions
    elif params.use_random:
        task_output("Randomizing transcripts")
        randomization_method = prepare_coordinates.RandomizeWithoutRegions
    else:
        task_output("No randomization selected, sorting input bed file")
        randomization_method = prepare_coordinates.NonRandom
    randomization_method.apply(params)

    if params.exclude_drivers:
        task_output("Excluding positively selected genes")
        drivers_method = prepare_coordinates.GeneExclusions
    else:
        task_output("Retaining positively selected genes")
        drivers_method = prepare_coordinates.GeneExclusionsDisabled
    drivers_method.apply(params)

    task_output("Building protein complement")
    prepare_coordinates.BuildProteinComplement.apply(params)

    if params.use_ssb192:
        task_output("Preparing CDS coords with SSB192 mutrate")
        ssb192_method = prepare_coordinates.UseSSB192
    else:
        task_output("Preparing CDS coords without SSB192 mutrate")
        ssb192_method = prepare_coordinates.NotSSB192
    ssb192_method.apply(params)

    task_output("Building intra epitope CDS file")
    prepare_coordinates.BuildIntraEpitopesCDS.apply(params)

    task_output("Obtaining fasta regions")
    obtain_fasta_regions.ObtainFastaRegions.apply(params)

    task_output(
        "Compiling list of transcript:regions to estimate number of sites"
    )
    obtain_fasta_regions.GetTranscriptRegionsForSites.apply(params)

    if params.use_ssb192:
        task_output(
            "Estimating all theoretically possible 192 substitutions for "
            "target and non-target regions"
        )
        analysis.ComputeSSB192TheoreticalSubs.apply(params)

        task_output(
            "Computing sum over all possible sites in target and non-target "
            "regions"
        )
        analysis.SumPossibleAcrossRegions.apply(params)

        task_output(
            "Processing VEP annotated file to estimated 192 rate parameters"
        )
        analysis.FixSimulated.apply(params)
        task_output("Applying column corrections for alternative format")
        analysis.ColumnCorrect.apply(params)

        task_output("Performing context corrections")
        analysis.ContextCorrection.apply(params)

        task_output("Flagging calculations")
        analysis.FlagComputations.apply(params)

        task_output("Computing triplet counts")
        analysis.TripletCounts.apply(params)

        task_output("Applying site corrections")
        analysis.SiteCorrections.apply(params)
    else:
        # TODO: See line 243 - 306
        raise ValueError("Implement SSB7")

    task_output("Intersecting by frequency")
    intersect.IntersectByFrequency.apply(params)


def parse_genome_args():
    parser = argparse.ArgumentParser(description="Genome reference")

    parser.add_argument(
        "--reference",
        "-r",
        dest="genome_ref",
        type=str,
        help="GRCh37 or GRCh38",
        required=True,
    )

    ref = check_genome(parser.parse_args())

    return ref


def download_genome():
    ref = parse_genome_args()
    startup_output()
    soprano_root_dir = pathlib.Path(__file__).parent
    installers_dir = soprano_root_dir.joinpath("bash_installers")
    downloader_path = installers_dir.joinpath("download_reference.sh")
    data_dir = soprano_root_dir.joinpath("data")

    assert downloader_path.exists(), downloader_path
    assert data_dir.exists(), data_dir

    subprocess.run(
        ["bash", downloader_path.as_posix(), ref, data_dir.as_posix()]
    )


if __name__ == "__main__":
    main()
