from SOPRANO import obtain_fasta_regions


def test__get_target_fasta_regions(step_3_defs):
    paths, transcripts = step_3_defs
    obtain_fasta_regions._get_target_fasta_regions(paths, transcripts)
