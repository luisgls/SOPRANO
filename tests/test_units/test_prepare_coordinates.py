import pathlib

import pytest

import SOPRANO.prepare_coordinates as prep_coords
from SOPRANO.objects import AnalysisPaths, AuxiliaryPaths, TranscriptPaths


def tab_line(*args):
    return "\t".join([str(arg) for arg in args]) + "\n"


# Snippet from TCGA-05-4396-01A-21D-1855-08.annotated
mock_input_content = [
    tab_line(
        "10_118957090_G/C",
        "10:118957090",
        "C",
        "ENSG00000186795",
        "ENST00000334549",
        "Transcript",
        "missense_variant",
        "91",
        "91",
        "31",
        "V/L",
        "Gtg/Ctg",
        "-",
        "IMPACT=MODERATE;STRAND=1;SYMBOL=KCNK18;SYMBOL_SOURCE=HGNC;HGNC_ID=19439",
    ),
    tab_line(
        "10_16528530_G/A",
        "10:16528530",
        "A",
        "ENSG00000165983",
        "ENST00000378000",
        "Transcript",
        "synonymous_variant",
        "858",
        "612",
        "204",
        "R",
        "cgG/cgA",
        "-",
        "IMPACT=LOW;STRAND=1;SYMBOL=PTER;SYMBOL_SOURCE=HGNC;HGNC_ID=9590",
    ),
    tab_line(
        "10_22021940_G/T",
        "10:22021940",
        "T",
        "ENSG00000078403",
        "ENST00000307729",
        "Transcript",
        "missense_variant",
        "2509",
        "2331",
        "777",
        "Q/H",
        "caG/caT",
        "-",
        "IMPACT=MODERATE;STRAND=1;SYMBOL=MLLT10;SYMBOL_SOURCE=HGNC;HGNC_ID=16063",
    ),
    tab_line(
        "10_43015604_C/A",
        "10:43015604",
        "A",
        "ENSG00000234420",
        "ENST00000452075",
        "Transcript",
        "non_coding_transcript_exon_variant",
        "1991",
        "-",
        "-",
        "-",
        "-",
        "-",
        "IMPACT=MODIFIER;STRAND=-1;SYMBOL=ZNF37BP;SYMBOL_SOURCE=HGNC;HGNC_ID=13103",
    ),
    tab_line(
        "10_51363129_C/A",
        "10:51363129",
        "A",
        "ENSG00000225784",
        "ENST00000404618",
        "Transcript",
        "non_coding_transcript_exon_variant",
        "943",
        "-",
        "-",
        "-",
        "-",
        "-",
        "IMPACT=MODIFIER;STRAND=-1;SYMBOL=RP11-592B15.4;SYMBOL_SOURCE=Clone_based_vega_gene",
    ),
]

# Snippet from TCGA-05-4396.Expressed.IEDBpeps.SB.epitope.bed
mock_bed_content = [
    tab_line("ENST00000000233", 115, 124),
    tab_line("ENST00000000233", 164, 177),
    tab_line("ENST00000000412", 113, 124),
    tab_line("ENST00000001008", 27, 36),
    tab_line("ENST00000001008", 189, 198),
]

# Snippet from ensemble_transcript.length
mock_transcript_content = [
    tab_line("ENST00000000233", 543),
    tab_line("ENST00000000412", 834),
    tab_line("ENST00000000442", 1272),
    tab_line("ENST00000001008", 1380),
    tab_line("ENST00000001146", 1539),
]

# ensemble_transcript_protein.length
mock_protein_transcript_content = [
    tab_line("ENST00000000233", 180),
    tab_line("ENST00000000412", 277),
    tab_line("ENST00000000442", 423),
    tab_line("ENST00000001008", 459),
    tab_line("ENST00000001146", 512),
]

# Fictitious target regions for randomization
mock_target_regions = [tab_line("ENST00000000233", 500, 1000)]

# Used to build mock aux files
mock_genes2exclude = [tab_line("ENST00000001008")]


def check_expected_content(
    expected_content: list, written_content_path: pathlib.Path
):
    assert written_content_path.exists()

    with open(written_content_path, "r") as f:
        written_content = f.readlines()

    assert len(expected_content) == len(written_content)

    for e, w in zip(expected_content, written_content):
        assert e.strip() == w.strip()


@pytest.fixture
def test_files(tmp_path):
    inputs_dir: pathlib.Path = tmp_path.joinpath("inputs")
    transcripts_dir: pathlib.Path = tmp_path.joinpath("transcripts")
    tmpdir: pathlib.Path = tmp_path.joinpath("tmp")

    inputs_dir.mkdir(parents=True)
    transcripts_dir.mkdir(parents=True)
    tmpdir.mkdir(parents=True)

    anno_path = inputs_dir.joinpath("input.anno")
    bed_path = inputs_dir.joinpath("input.bed")
    targets_path = inputs_dir.joinpath("targets.bed")
    genes2exclude_path = inputs_dir.joinpath("genes2exclude.txt")

    trans_path = transcripts_dir.joinpath("transcript_length.txt")
    trans_prot_path = transcripts_dir.joinpath(
        "transcript_length_protein.length"
    )

    paths = AnalysisPaths(
        "test_data", bed_path, tmpdir, target_regions_path=targets_path
    )
    transcripts = TranscriptPaths(trans_path, trans_prot_path)
    auxiliaries = AuxiliaryPaths(genes2exclude_path)

    for _input_path, _input_content in zip(
        (
            anno_path,
            bed_path,
            trans_path,
            trans_prot_path,
            targets_path,
            genes2exclude_path,
        ),
        (
            mock_input_content,
            mock_bed_content,
            mock_transcript_content,
            mock_protein_transcript_content,
            mock_target_regions,
            mock_genes2exclude,
        ),
    ):
        with open(_input_path, "w") as f:
            f.writelines(_input_content)

    return paths, transcripts, auxiliaries


@pytest.mark.dependency(name="_filter_transcript_file")
def test__filter_transcript_file(test_files):
    paths, transcripts, auxiliaries = test_files

    expected_content = [
        tab_line("ENST00000000233", 543),
        tab_line("ENST00000000412", 834),
        tab_line("ENST00000001008", 1380),
    ]

    prep_coords._filter_transcript_file(
        paths.bed_path,
        transcripts.transcript_length,
        paths.filtered_transcript,
    )

    check_expected_content(expected_content, paths.filtered_transcript)


@pytest.mark.dependency(
    name="filter_trans_files", depends=["_filter_transcript_file"]
)
def test_filter_transcript_files(test_files):
    paths, transcripts, auxiliaries = test_files

    expected_trans_content = [
        tab_line("ENST00000000233", 543),
        tab_line("ENST00000000412", 834),
        tab_line("ENST00000001008", 1380),
    ]

    expected_trans_protein_content = [
        tab_line("ENST00000000233", 180),
        tab_line("ENST00000000412", 277),
        tab_line("ENST00000001008", 459),
    ]

    prep_coords.filter_transcript_files(paths, transcripts)

    check_expected_content(expected_trans_content, paths.filtered_transcript)
    check_expected_content(
        expected_trans_protein_content, paths.filtered_protein_transcript
    )


@pytest.mark.dependency(name="_define_excl_regs")
def test__define_excluded_regions_for_randomization(test_files):
    paths, transcripts, auxiliaries = test_files

    expected_content = mock_bed_content + [
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000412", 0, 2),
        tab_line("ENST00000001008", 0, 2),
        tab_line("ENST00000001008", 0, 2),
    ]

    prep_coords._define_excluded_regions_for_randomization(paths)

    check_expected_content(expected_content, paths.exclusions)


@pytest.mark.dependency(depends=["_define_excl_regs", "filter_trans_files"])
def test__sort_excluded_regions_for_randomization(test_files):
    paths, transcripts, auxiliaries = test_files

    # Every item has an "ENST<...>  0   2" pair
    # So after sorting, we expect that
    expected_sorted_exclusions = [
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000233", 115, 124),
        tab_line("ENST00000000233", 164, 177),
        tab_line("ENST00000000412", 0, 2),
        tab_line("ENST00000000412", 113, 124),
        tab_line("ENST00000001008", 0, 2),
        tab_line("ENST00000001008", 0, 2),
        tab_line("ENST00000001008", 27, 36),
        tab_line("ENST00000001008", 189, 198),
    ]

    prep_coords.filter_transcript_files(paths, transcripts)
    prep_coords._define_excluded_regions_for_randomization(paths)
    prep_coords._sort_excluded_regions_for_randomization(paths, seed=1234)

    check_expected_content(expected_sorted_exclusions, paths.exclusions_sorted)

    ENST00000000412_min_max_vals = [
        [int(e.split("\t")[1]), int(e.split("\t")[2])]
        for e in expected_sorted_exclusions
        if e.startswith("ENST00000000412")
    ]

    with open(paths.exclusions_shuffled, "r") as f:
        shuffled_content = f.readlines()

    for line in shuffled_content:
        line = line.strip()
        chrom, start, stop, *other = line.split("\t")
        if chrom == "ENST00000000412":
            start = int(start)
            stop = int(stop)
            for input_start, input_stop in ENST00000000412_min_max_vals:
                overlap = (start <= input_stop) and (input_start <= stop)
                assert not overlap


@pytest.mark.dependency(depends=["filter_trans_files"])
def test__randomize_with_target_file(test_files):
    paths, transcripts, auxiliaries = test_files

    prep_coords.filter_transcript_files(paths, transcripts)

    with open(paths.filtered_transcript, "r") as f:
        expected_trans_content = f.readlines()
        expected_trans_content.append(tab_line("ENST00000000233", 543))

    with open(paths.filtered_protein_transcript, "r") as f:
        expected_trans_prot_content = f.readlines()
        expected_trans_prot_content.append(tab_line("ENST00000000233", 180))

    prep_coords._randomize_with_target_file(paths, transcripts, seed=1234)

    check_expected_content(expected_trans_content, paths.filtered_transcript)
    check_expected_content(
        expected_trans_prot_content, paths.filtered_protein_transcript
    )

    # TODO: Complete


def test__non_randomized(test_files):
    paths, transcripts, auxiliaries = test_files
    # expected_content = [
    #     tab_line("ENST00000000233", 115, 124),
    #     tab_line("ENST00000000233", 164, 177),
    #     tab_line("ENST00000000412", 113, 124),
    #     tab_line("ENST00000001008", 27, 36),
    #     tab_line("ENST00000001008", 189, 198),
    # ]
    prep_coords._non_randomized(paths)

    # TODO: See docstring in _non_randomized
    # check_expected_content(expected_content, paths.exclusions_shuffled)


def test__exclude_positively_selected_genes_disabled(test_files):
    paths, transcripts, auxiliaries = test_files

    # Dummy data written to shuffled exclusions file:
    # When positively selected genes are disabled, should just copy this file
    paths.exclusions_shuffled.write_text(tab_line("chr1", 123, 456))

    prep_coords._exclude_positively_selected_genes_disabled(paths)

    with open(paths.exclusions_shuffled, "r") as f:
        expected_content = f.readlines()

    check_expected_content(expected_content, paths.epitopes)


def test__exclude_positively_selected_genes(test_files):
    paths, transcripts, auxiliaries = test_files

    # Write dummy data for "shuffle file"
    paths.exclusions_shuffled.write_text(
        tab_line("ENST00000000233", 164, 177)
        + tab_line("ENST00000001008", 113, 124)
    )

    # In the exclusions aux file, we have ENST00000001008
    expected_content = [tab_line("ENST00000000233", 164, 177)]
    prep_coords._exclude_positively_selected_genes(paths, auxiliaries)

    check_expected_content(expected_content, paths.epitopes)


def test__build_protein_complement_tmp(test_files):
    paths, transcripts, auxiliaries = test_files

    # A minimal (unsorted) epitope file will look something like
    paths.epitopes.write_text(
        tab_line("ENST00000000233", 400, 500)
        + tab_line("ENST00000000233", 100, 150)
    )

    # A minimal filtered protein transcript file will look something like
    paths.filtered_protein_transcript.write_text(
        tab_line("ENST00000000233", 500)
        + tab_line("ENST00000000412", 277)
        + tab_line("ENST00000001008", 459)
    )

    # The expected complement file will then be
    expected_tmp_content = [
        tab_line("ENST00000000233", 0, 100),
        tab_line("ENST00000000233", 150, 400),
        tab_line("ENST00000000412", 0, 277),
        tab_line("ENST00000001008", 0, 459),
    ]

    # NOTE: the length in the filtered protein must be >= the upperbound of
    # the largest stop position in the epitopes file
    # (for the associated chroms), otherwise
    # bedtools will raise a warning which is undetected via the subprocess

    # NOTE: In the absence of start, stop in the epitopes file,
    # chromosomes appearing in the filtered file will automatically be
    # provided the complementary interval [0, length]

    prep_coords.get_protein_complement(paths)
    check_expected_content(expected_tmp_content, paths.intra_epitopes_prot_tmp)

    # The final step of the calculation essentially finds word regex matches
    # with the initial input. Therefore, in this case, the expected content are
    # the lines starting with ENST00000000233
    check_expected_content(expected_tmp_content[:2], paths.intra_epitopes_prot)
