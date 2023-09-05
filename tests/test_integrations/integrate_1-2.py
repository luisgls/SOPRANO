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


def test_initial_parse(tmp_path):
    namespace = Namespace(
        name=name,
        input=input_file,
        bed_file=bed_file,
        output=tmp_path,
        bed_regions=None,
        use_ssb192=True,
        exclude_drivers=False,
        transcript=objects.EnsemblTranscripts.transcript_length,
        protein_transcript=objects.EnsemblTranscripts.protein_transcript_length,
    )

    params = objects.Parameters.from_namespace(namespace)

    main(namespace)

    assert params.filtered_transcript.exists()
    assert params.filtered_protein_transcript.exists()
