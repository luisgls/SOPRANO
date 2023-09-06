import pathlib


def tab_line(*args):
    return "\t".join([str(arg) for arg in args]) + "\n"


def check_expected_content(
    expected_content: list, written_content_path: pathlib.Path
):
    assert written_content_path.exists()

    with open(written_content_path, "r") as f:
        written_content = f.readlines()

    assert len(expected_content) == len(written_content)

    for e, w in zip(expected_content, written_content):
        assert e.strip() == w.strip()
