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
genome_ref = "grch37"


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
        genome_ref=genome_ref,
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

    # Check that fasta files are built
    assert params.epitopes_cds_fasta.exists()
    assert params.intra_epitopes_cds_fasta.exists()

    assert params.epitopes_trans_regs.exists()
    assert params.intra_epitopes_trans_regs.exists()

    # Check theoretical site estimates computed
    assert params.epitopes_trans_regs.exists()
    assert params.intra_epitopes_trans_regs.exists()

    # Check that site summations have been performed
    assert params.epitopes_trans_regs_sum.exists()
    assert params.intra_epitopes_trans_regs_sum.exists()

    # Check col corrections have been applied to annotated input file
    assert params.sim_fixed.exists()
    assert params.col_corrected.exists()

    # Check contextualised calculation
    assert params.contextualised.exists()

    # Check flagged calculation
    assert params.flagged.exists()

    # Check triplet founts
    assert params.triplet_counts.exists()

    # Check that site corrections have been computed
    assert params.final_epitope_corrections.exists()
    assert params.final_intra_epitope_corrections.exists()

    # Check that intersect by frequency has been performed
    assert params.epitope_nans.exists()
    assert params.intra_epitope_nans.exists()

    # Check mutation counts
    assert params.variants_silent.exists()
    assert params.variants_nonsilent.exists()
    assert params.variants_missense.exists()
    assert params.variants_intronic.exists()
    assert params.raw_silent_count.exists()
    assert params.raw_nonsilent_count.exists()
    assert params.raw_missense_count.exists()
    assert params.in_silent_count.exists()
    assert params.in_nonsilent_count.exists()
    assert params.in_missense_count.exists()
    assert params.out_silent_count.exists()
    assert params.out_nonsilent_count.exists()
    assert params.out_missense_count.exists()

    # Check epitopes data file is build
    assert params.data_epitopes.exists()

    # Check intron rate has been computed
    assert params.intron_rate.exists()
