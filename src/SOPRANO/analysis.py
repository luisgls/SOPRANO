import pathlib

from SOPRANO.objects import AnalysisPaths, GenomePaths, Parameters
from SOPRANO.pipeline_utils import (
    MissingDataError,
    _PipelineComponent,
    is_empty,
)
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


def _build_flag_file(paths: AnalysisPaths):
    """
    Implements:

    paste $NAME.tmp $NAME.tmp.bed |
        cut -f6,14 - |  awk -F "/" '{FS="/"}{OFS="\t"}{print $1,$2}' |
            awk -F "" '{FS=""}{OFS="\t"}{if( ($1==$6) && ($3!="-") )
                {print "GOOD"}else{print "FAIL"}}' > $NAME.flag

    :param paths:
    """

    # Note: We had erroneous behaviour using the chain of "-F" via sub pipes
    # so implement alternatively:
    subprocess_pipes.pipe(
        [
            "paste",
            paths.col_corrected.as_posix(),
            paths.contextualised.as_posix(),
        ],
        ["cut", "-f6,14"],
        ["tr", "/", r"\t"],
        [
            "awk",
            r'{if($1 == substr($3,2,1) && $2 != "-")'
            r'{print "GOOD"}else{print "FAIL"}}',
        ],
        output_path=paths.flagged,
    )


class FlagComputations(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (params.col_corrected, params.contextualised)
        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        FlagComputations.check_ready(params)
        _build_flag_file(params)


def _initial_triplet_counts(paths: AnalysisPaths):
    r"""
    Implements:
    paste $NAME.tmp $NAME.tmp.bed $NAME.flag |
        awk '{if($15=="GOOD"){print $0}}' - | cut -f6,14 - |
            sort -k2,2 |uniq -c | sed 's/^ \+//g' | sort -k1,1 -n |
                sed 's/ /\t/g' | awk '{OFS="\t"}{print $3,$2,$1}' |
                    sed -e 's/\t[A-Z]\//_/g' > $NAME.finalVEP.triplets.counts
    :param paths:
    """

    subprocess_pipes.pipe(
        [
            "paste",
            paths.col_corrected.as_posix(),
            paths.contextualised.as_posix(),
            paths.flagged.as_posix(),
        ],
        ["awk", '{if($15=="GOOD"){print $0}}', "-"],
        ["cut", "-f6,14"],
        ["sort", "-k2,2"],
        ["uniq", "-c"],
        ["sed", r"s/^ \+//g"],
        ["sort", "-k1,1", "-n"],
        ["sed", r"s/ /\t/g"],
        ["awk", r'{OFS="\t"}{print $3,$2,$1}'],
        ["sed", "-e", r"s/\t[A-Z]\//_/g"],
        output_path=paths.triplet_counts,
    )


class TripletCounts(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (params.col_corrected, params.contextualised, params.flagged)
        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        TripletCounts.check_ready(params)
        _initial_triplet_counts(params)


class SOPRANOError(Exception):
    pass


def _check_triplet_counts(paths: AnalysisPaths):
    """
    Implements:

    muts=`wc -l $NAME.tmp.bed | awk '{ print $1 }'`
    back=`wc -l $NAME.finalVEP.triplets.counts | awk '{ print $1 }'`
    fails=`grep -c "FAIL" $NAME.flag`
    echo "Rate parameter file $NAME.finalVEP.triplets.counts has $back lines
    of data."
    rm $NAME.tmp $NAME.tmp.bed $NAME.flag
    echo "Proccesed $muts mutations from VEP file. $fails mutations were
    discarded (indels or reference conflicts)"
    if [ "$back" -lt 7 ]
    then
            echo "Very few rate parameters to analyse data, please provide
            rate parameter file with option -c"
            exit 1
        fi
    else
    echo "$NAME.finalVEP.triplets.counts is empty."
    exit 1

    :param paths:
    :return:
    """

    if is_empty(paths.triplet_counts):
        raise SOPRANOError(f"Triplet counts are empty: {paths.triplet_counts}")

    mutations = subprocess_pipes.pipe(
        ["wc", "-l", paths.flagged.as_posix()], ["awk", "{ print $1 }"]
    )
    back = subprocess_pipes.pipe(
        ["wc", "-l", paths.triplet_counts.as_posix()], ["awk", "{ print $1 }"]
    )
    fails = subprocess_pipes.pipe(
        ["grep", "-c", "FAIL", paths.flagged.as_posix()]
    )

    print(
        f"Rate parameter file {paths.triplet_counts.as_posix()} has {back} "
        f"lines of data.\n"
        f"Processed {mutations} from VEP file. {fails} mutations were "
        f"discarded (indels or reference conflicts)."
    )

    if int(back) < 7:
        raise SOPRANOError(
            "Too few rate parameters available to analyse data."
        )


def _correct_from_total_sites(paths: AnalysisPaths):
    """
    Implements:
    perl $BASEDIR/scripts/correct_update_epitope_sitesV3.pl
        $TMP/$NAME.listA.totalsites $NAME.finalVEP.triplets.counts >
            $TMP/$NAME.final_corrected_matrix_A.txt
    perl $BASEDIR/scripts/correct_update_epitope_sitesV3.pl
        $TMP/$NAME.listB.totalsites $NAME.finalVEP.triplets.counts >
            $TMP/$NAME.final_corrected_matrix_B.txt

    NOTE: listA.total sites is the trans_regs_sum for epitopes
    and listB.total sites is the trans_regs_sum for intra epitopes
    in the parameter defs.

    :param paths:
    """

    perl_script = SCRIPTS_DIR.joinpath("correct_update_epitope_sitesV3.pl")

    subprocess_pipes.pipe(
        [
            "perl",
            perl_script.as_posix(),
            paths.epitopes_trans_regs_sum.as_posix(),
            paths.triplet_counts.as_posix(),
        ],
        output_path=paths.final_epitope_corrections,
    )

    subprocess_pipes.pipe(
        [
            "perl",
            perl_script.as_posix(),
            paths.intra_epitopes_trans_regs_sum.as_posix(),
            paths.triplet_counts.as_posix(),
        ],
        output_path=paths.final_intra_epitope_corrections,
    )


class SiteCorrections(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.epitopes_trans_regs_sum,
            params.intra_epitopes_trans_regs_sum,
            params.triplet_counts,
        )

        for p in paths:
            if not p.exists():
                raise MissingDataError(p)

        _check_triplet_counts(params)

    @staticmethod
    def apply(params: Parameters):
        SiteCorrections.check_ready(params)
        _correct_from_total_sites(params)
