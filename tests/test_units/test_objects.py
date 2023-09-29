import pathlib

from SOPRANO.core import objects


def test_Transcripts():
    defaults = objects.TranscriptPaths.defaults()

    assert defaults.transcript_length.exists()
    assert defaults.protein_transcript_length.exists()
    assert defaults.transcript_fasta.exists()


def test_GRCh37():
    assert objects.GenomePaths.GRCh37().sizes.exists()


def test_GRCh38():
    assert objects.GenomePaths.GRCh38().sizes.exists()


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
