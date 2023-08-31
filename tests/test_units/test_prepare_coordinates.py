import pathlib
import tempfile

import pytest

import SOPRANO.prepare_coordinates as prep_coords
from SOPRANO.objects import AnalysisPaths, TranscriptPaths


def tab_line(*args):
    return "\t".join([str(arg) for arg in args]) + "\n"


@pytest.mark.dependency(name="_filter_transcript_file")
def test__filter_transcript_file():
    mock_bed_content = [
        "ENST00000000003 9      180\n",
        "ENST00000000001 9       11\n",
        "ENST00000000002 99      18\n",
    ]

    mock_transcript_content = [
        "ENST00000000001 2738\n",
        "ENST00000000002 12\n",
    ]

    expected_content = ["ENST00000000001", "ENST00000000002"]

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)

        bed_path = tmp_dir.joinpath("file.bed")
        trans_path = tmp_dir.joinpath("trans.length")
        filt_path = tmp_dir.joinpath("trans_filt.length")

        with open(bed_path, "w") as b:
            b.writelines(mock_bed_content)

        with open(trans_path, "w") as t:
            t.writelines(mock_transcript_content)

        prep_coords._filter_transcript_file(bed_path, trans_path, filt_path)

        assert filt_path.exists(), filt_path

        with open(filt_path, "r") as f:
            written_content = f.readlines()

        for e, w in zip(expected_content, written_content):
            assert e == w.strip()


@pytest.mark.dependency(depends=["_filter_transcript_file"])
def test_filter_transcript_files():
    mock_bed_content = [
        "ENST00000000003 9      180\n",
        "ENST00000000001 9       11\n",
        "ENST00000000002 99      18\n",
    ]

    mock_transcript_content = [
        "ENST00000000001 2738\n",
        "ENST00000000002 12\n",
    ]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)

        bed_path = tmpdir.joinpath("bed")

        trans_path = tmpdir.joinpath("trans.length")
        protein_path = tmpdir.joinpath("protein.length")

        for path, content in zip(
            (bed_path, trans_path, protein_path),
            (
                mock_bed_content,
                mock_transcript_content,
                mock_transcript_content,
            ),
        ):
            with open(path, "w") as f:
                f.writelines(content)

        paths = AnalysisPaths("test", bed_path, tmpdir)
        transcripts = TranscriptPaths(trans_path, protein_path)

        prep_coords.filter_transcript_files(paths, transcripts)

        filt_trans_path = paths.filtered_transcript
        filt_protein_path = paths.filtered_protein_transcript
        assert filt_trans_path.exists()
        assert filt_protein_path.exists()


def test__define_excluded_regions_for_randomization():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]
    expected_bed_content = [
        tab_line("chr1", 800, 1000),
        tab_line("chr1", 80, 180),
        tab_line("chr1", 1, 10),
        tab_line("chr1", 750, 10000),
        tab_line("chr1", 0, 2),
        tab_line("chr1", 0, 2),
        tab_line("chr1", 0, 2),
        tab_line("chr1", 0, 2),
    ]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)
        name = "test"
        bed_path = tmpdir.joinpath("bed.file")

        with open(bed_path, "w") as f:
            f.writelines(mock_bed_content)

        result_path = tmpdir.joinpath(f"{name}.exclusion.ori")

        paths = AnalysisPaths(name, bed_path, tmpdir)

        prep_coords._define_excluded_regions_for_randomization(paths)

        assert result_path.exists()

        with open(result_path, "r") as f:
            written_content = f.readlines()

        assert len(expected_bed_content) == len(written_content)

        for e, w in zip(expected_bed_content, written_content):
            assert e.strip() == w.strip()


def test__sort_excluded_regions_for_randomization():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]

    mock_filtered_protein = [
        tab_line("chr1", 10000),
        tab_line("chr2", 8000),
        tab_line("chr3", 5000),
        tab_line("chr2", 2000),
    ]

    input_exclusions = [
        tab_line("chr1", 300, 400),
        tab_line("chr1", 200, 250),
    ]

    # Default bedsort is by chrom then start pos
    expected_sorted_exclusions = [
        tab_line("chr1", 200, 250),
        tab_line("chr1", 300, 400),
    ]
    chr1_exclusions = [[200, 250], [300, 400]]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)
        name = "test"
        bed_path = tmpdir.joinpath("bed.file")
        paths = AnalysisPaths(name, bed_path, tmpdir)

        with open(paths.bed_path, "w") as f_bed:
            f_bed.writelines(mock_bed_content)

        with open(paths.exclusions, "w") as f_excl:
            f_excl.writelines(input_exclusions)

        with open(paths.filtered_protein_transcript, "w") as f_prot:
            f_prot.writelines(mock_filtered_protein)

        prep_coords._sort_excluded_regions_for_randomization(paths, seed=1234)

        assert paths.exclusions_sorted.exists()

        with open(paths.exclusions_sorted, "r") as f_sort:
            written_sorted_exclusions = f_sort.readlines()

        for e, w in zip(expected_sorted_exclusions, written_sorted_exclusions):
            assert e.strip() == w.strip()

        with open(paths.exclusions_shuffled, "r") as f_shuf:
            written_shuffled = f_shuf.readlines()

        for line in written_shuffled:
            line = line.strip()
            chrom, start, stop, *other = line.split("\t")
            if chrom == "chr1":
                start = int(start)
                stop = int(stop)
                for input_start, input_stop in chr1_exclusions:
                    overlap = (start <= input_stop) and (input_start <= stop)
                    assert not overlap


def test__randomize_with_target_file():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]

    mock_filtered_transcript = [
        tab_line("chr1", 10000),
        tab_line("chr2", 8000),
        tab_line("chr3", 5000),
        tab_line("chr2", 2000),
    ]

    mock_transcript_content = [
        tab_line("chr1", 2738),
        tab_line("chr2", 12),
    ]

    targets = [tab_line("chr1", 500, 1000), tab_line("chr1", 2000, 2500)]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)

        paths = AnalysisPaths(
            "test",
            tmpdir.joinpath("test.bed"),
            tmpdir,
            target_regions_path=tmpdir.joinpath("targets.bed"),
        )

        transcripts = TranscriptPaths(
            tmpdir.joinpath("trans.length"),
            tmpdir.joinpath("trans_prot.length"),
        )

        for path, lines in zip(
            (
                paths.bed_path,
                paths.filtered_transcript,
                paths.filtered_protein_transcript,
                paths.target_regions_path,
                transcripts.transcript_length,
                transcripts.protein_transcript_length,
            ),
            (
                mock_bed_content,
                mock_filtered_transcript,
                mock_filtered_transcript,
                targets,
                mock_transcript_content,
                mock_transcript_content,
            ),
        ):
            with open(path, "w") as f:
                f.writelines(lines)

        prep_coords._randomize_with_target_file(paths, transcripts, seed=1234)

        with open(paths.filtered_protein_transcript, "r") as f:
            written_filtered_lines = f.readlines()

        expected_filtered_lines = mock_filtered_transcript + [
            tab_line("chr1", 2738)
        ]

        for e, w in zip(expected_filtered_lines, written_filtered_lines):
            assert e.strip() == w.strip()


def test__non_randomized():
    mock_bed_content = [
        tab_line("chr1", 800, 1000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
    ]

    expected_content = [
        tab_line("chr1", 1, 10, 24),
        tab_line("chr1", 750, 10000, 24),
        tab_line("chr1", 80, 180, 24),
        tab_line("chr1", 800, 1000, 24),
    ]

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)

        paths = AnalysisPaths("test", tmpdir.joinpath("test.bed"), tmpdir)

        with open(paths.bed_path, "w") as f:
            f.writelines(mock_bed_content)

        prep_coords._non_randomized(paths)

        assert paths.exclusions_shuffled.exists()

        with open(paths.exclusions_shuffled, "r") as f:
            written = f.readlines()

        for e, w in zip(expected_content, written):
            assert e.strip() == w.strip()
