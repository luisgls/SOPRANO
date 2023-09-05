import pathlib

import pytest

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
        "test_data", bed_path, tmpdir, target_regions=targets_path
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


@pytest.fixture
def minimal_epitopes(tmp_path):
    epi_dir: pathlib.Path = tmp_path.joinpath("epitopes")
    epi_dir.mkdir(parents=True)

    paths = AnalysisPaths(
        "test_epitopes", tmp_path.joinpath("mock.bed"), epi_dir
    )

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

    return paths
