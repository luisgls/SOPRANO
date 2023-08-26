import pathlib
import tempfile

import SOPRANO.prepare_coordinates as prep_coords


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

    # expected_content = ["ENST00000000001\n", "ENST00000000002\n"]

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

        # TODO: Check output from filtered lines is as
        #       expected


test__filter_transcript_file()
