import pathlib

from SOPRANO.objects import Parameters
from SOPRANO.pipeline_utils import MissingDataError, _PipelineComponent
from SOPRANO.sh_utils import subprocess_pipes

SOPRANO_ROOT = pathlib.Path(__file__).parent
SCRIPTS_DIR = SOPRANO_ROOT.joinpath("scripts")


def _compute_theoretical_subs(
    cds_fasta: pathlib.Path, trans_regs: pathlib.Path
):
    """

    Implement

    $BASEDIR/scripts/calculate_sites_signaturesLZ_192.pl
        $TMP/$NAME.epitopes_cds.fasta $TMP/$NAME.listA > $TMP/$NAME.listA.sites

    Estimates all theoretical possible 192 subsitutions in target and
    non-target regions

    :param cds_fasta: path to epitopes_cds of intra_epitopes_cds fasta file
    :param trans_regs: list of transcript:regions
    """

    perl_path = SCRIPTS_DIR.joinpath("calculate_sites_signaturesLZ_192.pl")

    subprocess_pipes.pipe(
        [perl_path.as_posix(), cds_fasta.as_posix(), trans_regs],
        output_path=trans_regs,
        overwrite=True,
    )


class ComputeSSB192TheoreticalSubs(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.epitopes_cds_fasta,
            params.intra_epitopes_cds_fasta,
            params.epitopes_trans_regs,
            params.intra_epitopes_trans_regs,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ComputeSSB192TheoreticalSubs.check_ready(params)
        _compute_theoretical_subs(
            params.epitopes_cds_fasta, params.epitopes_trans_regs
        )
        _compute_theoretical_subs(
            params.intra_epitopes_cds_fasta, params.intra_epitopes_trans_regs
        )


def _sum_possible_across_region(
    trans_regs: pathlib.Path, sum_trans_regs: pathlib.Path
):
    """
    Implements:
    awk '{print "test_"$2"_"$3"\t0\t1\t"$0}' $TMP/$NAME.listA.sites |
        sortBed -i stdin | mergeBed -i stdin -c 7,8 -o sum,sum |
            cut -f1,4,5 > $TMP/$NAME.listA.totalsites

    :param trans_regs: list of transcript:regions
    :param sum_trans_regs: output path for summation
    """

    subprocess_pipes.pipe(
        ["awk", '{print "test_"$2"_"$3"\t0\t1\t"$0}', trans_regs.as_posix()],
        ["sortBed", "-i", "stdin"],
        ["mergeBed", "-i", "stdin", "-c", "7,8", "-o", "sum,sum"],
        ["cut", "-f1,4,5"],
        output_path=sum_trans_regs,
    )


class SumPossibleAcrossRegions(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (params.epitopes_trans_regs, params.intra_epitopes_trans_regs)
        for p in paths:
            if not p.exists():
                raise MissingDataError(p)

    @staticmethod
    def apply(params: Parameters):
        SumPossibleAcrossRegions.check_ready(params)
        _sum_possible_across_region(
            params.epitopes_trans_regs, params.epitopes_trans_regs_sum
        )
        _sum_possible_across_region(
            params.intra_epitopes_trans_regs,
            params.intra_epitopes_trans_regs_sum,
        )
