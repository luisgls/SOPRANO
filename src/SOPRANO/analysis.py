import pathlib

from SOPRANO.objects import AnalysisPaths, GenomePaths, Parameters
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


def _fix_simulated(paths: AnalysisPaths):
    """
    Implements:

    perl $BASEDIR/scripts/fixsimulated.pl $FILE > $NAME

    :param paths:
    """

    perl_script = SCRIPTS_DIR.joinpath("fixsimulated.pl")

    subprocess_pipes.pipe(
        ["perl", perl_script, paths.input_path], output_path=paths.sim_fixed
    )


class FixSimulated(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        FixSimulated.check_ready(params)
        _fix_simulated(params)


def _col_correct(paths: AnalysisPaths):
    """
    Implements:

    cut -f1,5,7,11,12 $NAME | grep -v "#" | sed 's/_/\t/1' | sed 's/_/\t/1' |
        awk -F"\t"
            '{OFS="\t"}{if($6!="-"&&length($3)==3){print $1,$2-1,$2,$0}}' >
                $NAME.tmp

    :param paths:
    """

    # NOTE: Removed -F"\t" component - seems unneeded, and breaks the
    # subprocess pipe
    subprocess_pipes.pipe(
        ["cut", "-f1,5,7,11,12", paths.sim_fixed.as_posix()],
        ["grep", "-v", "#"],
        ["sed", "s/_/\t/1"],
        ["sed", "s/_/\t/1"],
        [
            "awk",
            '{OFS="\t"}{if($6!="-" && length($3)==3){print $1,$2-1,$2,$0}}',
        ],
        output_path=paths.col_corrected,
    )


class ColumnCorrect(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        if not params.sim_fixed.exists():
            raise MissingDataError(params.sim_fixed)

    @staticmethod
    def apply(params: Parameters):
        ColumnCorrect.check_ready(params)
        _col_correct(params)


def _context_correction(paths: AnalysisPaths, genomes: GenomePaths):
    """
    Implements:
    bedtools slop -i $NAME.tmp -b 1 -g $GENOME |
        bedtools getfasta -fi $FASTA -bed stdin -tab |
            sed 's/:/\t/g' | sed 's/-/\t/g' > $NAME.tmp.bed

    :param paths:
    :return:
    """

    subprocess_pipes.pipe(
        [
            "bedtools",
            "slop",
            "-i",
            paths.col_corrected.as_posix(),
            "-b",
            "1",
            "-g",
            genomes.sizes.as_posix(),
        ],
        [
            "bedtools",
            "getfasta",
            "-fi",
            genomes.fasta.as_posix(),
            "-bed",
            "stdin",
            "-tab",
        ],
        ["sed", "s/:/\t/g"],
        ["sed", "s/-/\t/g"],
        output_path=paths.contextualised,
    )


class ContextCorrection(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.col_corrected,
            params.genomes.fasta,
            params.genomes.sizes,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ContextCorrection.check_ready(params)
        _context_correction(params, params.genomes)
