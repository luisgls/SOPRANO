import pathlib
from argparse import Namespace

import SOPRANO
from SOPRANO import objects
from SOPRANO.run_local_ssb_selection import main

SOPRANO_ROOT = pathlib.Path(SOPRANO.__file__).parent
DATA_DIR = SOPRANO_ROOT.joinpath("data")
BIO_DIR = SOPRANO_ROOT.joinpath("immunopeptidomes").joinpath("human")
EXAMPLES_DIR = SOPRANO_ROOT.joinpath("examples")

input_file = EXAMPLES_DIR.joinpath("TCGA-05-4396-01A-21D-1855-08.annotated")
bed_file = BIO_DIR.joinpath("TCGA-05-4396.Expressed.IEDBpeps.SB.epitope.bed")
name = "TCGA-05-4396"


def test_pipeline(tmp_path):
    namespace = Namespace(
        analysis_name=name,
        input_path=input_file,
        bed_path=bed_file,
        cache_dir=tmp_path,
        target_regions=None,
        use_ssb192=True,
        use_random=False,
        exclude_drivers=False,
        seed=-1,
        transcript=objects.EnsemblTranscripts.transcript_length,
        protein_transcript=objects.EnsemblTranscripts.protein_transcript_length,
        transcript_ids=objects.EnsemblTranscripts.transcript_fasta,
    )

    params = objects.Parameters.from_namespace(namespace)

    main(namespace)

    # Check filtered transcripts have been built
    assert params.filtered_transcript.exists()
    assert params.filtered_protein_transcript.exists()

    # Check drivers have (not) been excluded
    assert params.epitopes.exists()

    # Check epitope files been produced
    assert params.epitopes.exists()
    assert params.epitopes_cds.exists()

    # Check complement files been produced
    assert params.intra_epitopes_prot.exists()
    assert params.intra_epitopes_cds.exists()
    assert params.intra_epitopes_prot.exists()
