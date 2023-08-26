import pathlib
import tempfile

import pytest

import SOPRANO.prepare_coordinates as prep_coords


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

        prep_coords._filter_transcript_file(bed_path, trans_path, tmp_dir)

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

        prep_coords.filter_transcript_files(
            bed_path, trans_path, protein_path, tmpdir
        )

        filt_trans_path = tmpdir.joinpath("trans_filt.length")
        filt_protein_path = tmpdir.joinpath("protein_filt.length")

        assert filt_trans_path.exists()
        assert filt_protein_path.exists()
