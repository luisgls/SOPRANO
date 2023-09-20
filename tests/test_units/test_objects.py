import pathlib

from SOPRANO import objects


def test__data_dir():
    assert objects._data_dir().exists()


def test_Transcripts():
    assert objects.EnsemblTranscripts.transcript_length.exists()
    assert objects.EnsemblTranscripts.protein_transcript_length.exists()


def test_GRCh37():
    assert objects.GRCh37.sizes.exists()
    # assert objects.GRCh37.fasta.exists()


def test_cache_path_builder():
    test_dir = pathlib.Path.cwd()

    name = "foo"
    extensions = ("bar", "spam", "eggs")

    no_exs = test_dir.joinpath(name)
    with_exs = objects.cache_path_builder(test_dir, name, *extensions)

    assert no_exs.parent == with_exs.parent
    assert no_exs.name != with_exs

    for _ in range(len(extensions)):
        with_exs = with_exs.with_suffix("")

    assert no_exs == with_exs
