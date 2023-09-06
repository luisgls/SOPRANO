import pytest
from _test_utils import tab_line

from SOPRANO.objects import AnalysisPaths, TranscriptPaths

mock_tcga_epitopes_content = [
    tab_line("ENST00000000233", 115, 124),
    tab_line("ENST00000000233", 164, 177),
    tab_line("ENST00000000233", 37, 46),
    tab_line("ENST00000000233", 9, 18),
    tab_line("ENST00000000233", 98, 111),
    tab_line("ENST00000000412", 113, 128),
    tab_line("ENST00000001008", 189, 198),
    tab_line("ENST00000001008", 258, 267),
    tab_line("ENST00000001008", 27, 36),
    tab_line("ENST00000001008", 362, 380),
]

mock_tcga_epitopes_cds_content = [
    tab_line("ENST00000000233", 339, 375, "ENST00000000233", 115, 124),
    tab_line("ENST00000000233", 486, 534, "ENST00000000233", 164, 177),
    tab_line("ENST00000000233", 105, 141, "ENST00000000233", 37, 46),
    tab_line("ENST00000000233", 21, 57, "ENST00000000233", 9, 18),
    tab_line("ENST00000000233", 288, 336, "ENST00000000233", 98, 111),
    tab_line("ENST00000000412", 333, 387, "ENST00000000412", 113, 128),
    tab_line("ENST00000001008", 561, 597, "ENST00000001008", 189, 198),
    tab_line("ENST00000001008", 768, 804, "ENST00000001008", 258, 267),
    tab_line("ENST00000001008", 75, 111, "ENST00000001008", 27, 36),
    tab_line("ENST00000001008", 1080, 1143, "ENST00000001008", 362, 380),
]

mock_tcga_intra_epitopes_cds_content = [
    tab_line("ENST00000000233", 0, 21),
    tab_line("ENST00000000233", 57, 105),
    tab_line("ENST00000000233", 141, 288),
    tab_line("ENST00000000233", 336, 339),
    tab_line("ENST00000000233", 375, 486),
    tab_line("ENST00000000233", 534, 543),
    tab_line("ENST00000000412", 0, 333),
    tab_line("ENST00000000412", 387, 834),
    tab_line("ENST00000001008", 0, 75),
    tab_line("ENST00000001008", 111, 561),
]

mock_tcga_inta_epitopes_prot_content = [
    tab_line("ENST00000000233", 0, 9),
    tab_line("ENST00000000233", 18, 37),
    tab_line("ENST00000000233", 46, 98),
    tab_line("ENST00000000233", 111, 115),
    tab_line("ENST00000000233", 124, 164),
    tab_line("ENST00000000233", 177, 180),
    tab_line("ENST00000000412", 0, 113),
    tab_line("ENST00000000412", 128, 277),
    tab_line("ENST00000001008", 0, 27),
    tab_line("ENST00000001008", 36, 189),
]

mock_length_filt_content = [
    tab_line("ENST00000000233", 543),
    tab_line("ENST00000000412", 834),
    tab_line("ENST00000001008", 1380),
    tab_line("ENST00000001146", 1539),
    tab_line("ENST00000002125", 1326),
    tab_line("ENST00000002165", 1404),
    tab_line("ENST00000002829", 2358),
    tab_line("ENST00000003084", 4443),
    tab_line("ENST00000003100", 1530),
    tab_line("ENST00000003302", 3234),
]

mock_protein_length_filt_content = [
    tab_line("ENST00000000233", 180),
    tab_line("ENST00000000412", 277),
    tab_line("ENST00000001008", 459),
    tab_line("ENST00000001146", 512),
    tab_line("ENST00000002125", 441),
    tab_line("ENST00000002165", 467),
    tab_line("ENST00000002829", 785),
    tab_line("ENST00000003084", 1480),
    tab_line("ENST00000003100", 509),
    tab_line("ENST00000003302", 1077),
]

mock_transcript_id_content = [
    tab_line(">ENST00000448914"),
    tab_line("ACTGGGGGATACG"),
    tab_line(">ENST00000604642"),
    tab_line("GTGGATATAGTGTCTACGATTAC"),
    tab_line(">ENST00000603326"),
    tab_line("NNTGACTATGGTGCTAACTAC"),
    tab_line(">ENST00000604950"),
    tab_line("NNGTATTATGATTTTTGGACTGGTTATTATACC"),
    tab_line(">ENST00000603077"),
    tab_line("NNAGAATATTGTAATAGTACTACTTTCTATGCC"),
    tab_line(">ENST00000605284"),
    tab_line("GGTATAACTGGAACAAC"),
    tab_line(">ENST00000604446"),
    tab_line("GTGGATATAGTGTCTACGATTAC"),
    tab_line(">ENST00000603693"),
    tab_line("NNTGACTATGGTGCTAACTAC"),
]


@pytest.fixture
def step_3_defs(tmp_path):
    analysis_paths = AnalysisPaths(
        "step3", tmp_path.joinpath("foo"), tmp_path.joinpath("bar"), tmp_path
    )

    epitopes_path = analysis_paths.epitopes
    epitopes_cds_path = analysis_paths.epitopes_cds

    intra_path = analysis_paths.intra_epitopes
    intra_prot_path = analysis_paths.intra_epitopes_prot

    filt_path = analysis_paths.filtered_transcript
    filt_prot_path = analysis_paths.filtered_protein_transcript

    transcripts = TranscriptPaths(
        tmp_path.joinpath("spam"),
        tmp_path.joinpath("eggs"),
        tmp_path.joinpath("transcriptID.fasta"),
    )

    for path, content in zip(
        (
            epitopes_path,
            epitopes_cds_path,
            intra_path,
            intra_prot_path,
            filt_path,
            filt_prot_path,
            transcripts.transcript_fasta,
        ),
        (
            mock_tcga_epitopes_content,
            mock_tcga_epitopes_cds_content,
            mock_tcga_intra_epitopes_cds_content,
            mock_tcga_inta_epitopes_prot_content,
            mock_length_filt_content,
            mock_protein_length_filt_content,
            mock_transcript_id_content,
        ),
    ):
        with open(path, "w") as f:
            f.writelines(content)

    return analysis_paths, transcripts
