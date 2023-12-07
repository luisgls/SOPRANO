import pathlib

from SOPRANO.core.objects import AnalysisPaths, GenomePaths
from SOPRANO.utils.path_utils import Directories, is_empty
from SOPRANO.utils.sh_utils import pipe


def _compute_theoretical_subs_192(
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

    perl_path = Directories.scripts("calculate_sites_signaturesLZ_192.pl")

    pipe(
        [
            "perl",
            perl_path.as_posix(),
            cds_fasta.as_posix(),
            trans_regs.as_posix(),
        ],
        output_path=trans_regs,
        overwrite=True,
    )


def _compute_theoretical_subs_7(
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

    perl_path = Directories.scripts("calculate_sites_signaturesLZ.pl")

    pipe(
        [
            "perl",
            perl_path.as_posix(),
            cds_fasta.as_posix(),
            trans_regs.as_posix(),
        ],
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

    pipe(  # TODO: This appears to be a memory bottleneck ...
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

    perl_script = Directories.scripts("fixsimulated.pl")

    pipe(["perl", perl_script, paths.input_path], output_path=paths.sim_fixed)


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
    pipe(
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


def _context_correction(paths: AnalysisPaths, genomes: GenomePaths):
    """
    Implements:
    bedtools slop -i $NAME.tmp -b 1 -g $GENOME |
        bedtools getfasta -fi $FASTA -bed stdin -tab |
            sed 's/:/\t/g' | sed 's/-/\t/g' > $NAME.tmp.bed

    :param paths:
    :return:
    """

    pipe(
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
    pipe(
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

    pipe(
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


class SOPRANOError(Exception):
    pass


def _transform_192_to_7(paths: AnalysisPaths):
    r"""
    Implements:

    perl $BASEDIR/scripts/transform192to7.pl $NAME.finalVEP.triplets.counts
        $BASEDIR/data/final_translate_SSB192toSSB7 |
            awk -F "\t" '{OFS="\t"}{print $3,1,2,$2}' | sortBed -i stdin |
                mergeBed -i stdin -c 4 -o sum |
                    awk '{OFS="\t"}{print "Estimated",$1,$4}' |
                        sed 's/_/\//g' > tmp_to_7

    cp $NAME.finalVEP.triplets.counts $NAME.finalVEP.triplets192.counts
    mv tmp_to_7 $NAME.finalVEP.triplets.counts
    """

    perl_path = Directories.scripts("transform192to7.pl")
    translation_data_path = Directories.soprano_aux_files(
        "final_translate_SSB192toSSB7"
    )

    if is_empty(paths.triplet_counts):
        raise Warning(
            f"{paths.triplet_counts}is empty before attempting conversion"
        )

    assert perl_path.exists()
    assert translation_data_path.exists()

    import shutil

    shutil.copy(paths.triplet_counts, paths.triplet_counts.as_posix() + ".bkp")

    pipe(
        [
            "perl",
            perl_path.as_posix(),
            paths.triplet_counts.as_posix(),
            translation_data_path.as_posix(),
        ],
        ["awk", r'{OFS="\t"}{print $3,1,2,$2}'],  # removed -F "\t" bit
        ["sortBed", "-i", "stdin"],
        ["mergeBed", "-i", "stdin", "-c", "4", "-o", "sum"],
        ["awk", r'{OFS="\t"}{print "Estimated",$1,$4}'],
        ["sed", r"s/_/\//g"],
        output_path=paths.triplet_counts,
        overwrite=True,
    )


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

    mutations = pipe(
        ["wc", "-l", paths.flagged.as_posix()], ["awk", "{ print $1 }"]
    )
    back = pipe(
        ["wc", "-l", paths.triplet_counts.as_posix()], ["awk", "{ print $1 }"]
    )
    try:
        fails = pipe(["grep", "-c", "FAIL", paths.flagged.as_posix()])
    except RuntimeError:
        # If no FAIL are found, exit code is non-zero (confusingly!) so except
        fails = "0"

    print(
        f"Rate parameter file {paths.triplet_counts.as_posix()} has {back} "
        f"lines of data.\n"
        f"Processed {mutations} from annotated file. {fails} mutations were "
        f"discarded (indels or reference conflicts)."
    )

    if int(back) + 1 < 7:  # +1 due to lack of trailing \n in pipe output
        raise SOPRANOError(
            "Too few rate parameters available to analyse data."
        )


def _correct_from_total_sites_ssb192(paths: AnalysisPaths):
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

    perl_script = Directories.scripts("correct_update_epitope_sitesV3.pl")

    pipe(
        [
            "perl",
            perl_script.as_posix(),
            paths.epitopes_trans_regs_sum.as_posix(),
            paths.triplet_counts.as_posix(),
        ],
        output_path=paths.final_epitope_corrections,
    )

    pipe(
        [
            "perl",
            perl_script.as_posix(),
            paths.intra_epitopes_trans_regs_sum.as_posix(),
            paths.triplet_counts.as_posix(),
        ],
        output_path=paths.final_intra_epitope_corrections,
    )


def _correct_from_total_sites_ssb7(paths: AnalysisPaths):
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

    perl_script = Directories.scripts("correct_update_epitope_sitesV2.pl")

    pipe(
        [
            "perl",
            perl_script.as_posix(),
            paths.epitopes_trans_regs_sum.as_posix(),
            paths.triplet_counts.as_posix(),
        ],
        output_path=paths.final_epitope_corrections,
    )

    pipe(
        [
            "perl",
            perl_script.as_posix(),
            paths.intra_epitopes_trans_regs_sum.as_posix(),
            paths.triplet_counts.as_posix(),
        ],
        output_path=paths.final_intra_epitope_corrections,
    )
