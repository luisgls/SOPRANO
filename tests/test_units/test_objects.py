from SOPRANO.core import objects


def test_Transcripts():
    defaults = objects.TranscriptPaths.defaults()

    assert defaults.transcript_length.exists()
    assert defaults.protein_transcript_length.exists()


def test_GRCh37():
    assert objects.GenomePaths.GRCh37().sizes.exists()


def test_GRCh38():
    assert objects.GenomePaths.GRCh38().sizes.exists()
