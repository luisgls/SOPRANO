import pytest
from conftest import check_expected_content, mock_bed_content, tab_line

import SOPRANO.prepare_coordinates as prep_coords


@pytest.mark.dependency(name="_filter_transcript_file")
def test__filter_transcript_file(test_files):
    paths, transcripts, auxiliaries = test_files

    expected_content = [
        tab_line("ENST00000000233", 543),
        tab_line("ENST00000000412", 834),
        tab_line("ENST00000001008", 1380),
    ]

    prep_coords._filter_transcript_file(
        paths.bed_path,
        transcripts.transcript_length,
        paths.filtered_transcript,
    )

    check_expected_content(expected_content, paths.filtered_transcript)


@pytest.mark.dependency(
    name="filter_trans_files", depends=["_filter_transcript_file"]
)
def test_filter_transcript_files(test_files):
    paths, transcripts, auxiliaries = test_files

    expected_trans_content = [
        tab_line("ENST00000000233", 543),
        tab_line("ENST00000000412", 834),
        tab_line("ENST00000001008", 1380),
    ]

    expected_trans_protein_content = [
        tab_line("ENST00000000233", 180),
        tab_line("ENST00000000412", 277),
        tab_line("ENST00000001008", 459),
    ]

    prep_coords.filter_transcript_files(paths, transcripts)

    check_expected_content(expected_trans_content, paths.filtered_transcript)
    check_expected_content(
        expected_trans_protein_content, paths.filtered_protein_transcript
    )


@pytest.mark.dependency(name="_define_excl_regs")
def test__define_excluded_regions_for_randomization(test_files):
    paths, transcripts, auxiliaries = test_files

    expected_content = mock_bed_content + [
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000412", 0, 2),
        tab_line("ENST00000001008", 0, 2),
        tab_line("ENST00000001008", 0, 2),
    ]

    prep_coords._define_excluded_regions_for_randomization(paths)

    check_expected_content(expected_content, paths.exclusions)


@pytest.mark.dependency(depends=["_define_excl_regs", "filter_trans_files"])
def test__sort_excluded_regions_for_randomization(test_files):
    paths, transcripts, auxiliaries = test_files

    # Every item has an "ENST<...>  0   2" pair
    # So after sorting, we expect that
    expected_sorted_exclusions = [
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000233", 0, 2),
        tab_line("ENST00000000233", 115, 124),
        tab_line("ENST00000000233", 164, 177),
        tab_line("ENST00000000412", 0, 2),
        tab_line("ENST00000000412", 113, 124),
        tab_line("ENST00000001008", 0, 2),
        tab_line("ENST00000001008", 0, 2),
        tab_line("ENST00000001008", 27, 36),
        tab_line("ENST00000001008", 189, 198),
    ]

    prep_coords.filter_transcript_files(paths, transcripts)
    prep_coords._define_excluded_regions_for_randomization(paths)
    prep_coords._sort_excluded_regions_for_randomization(paths, seed=1234)

    check_expected_content(expected_sorted_exclusions, paths.exclusions_sorted)

    ENST00000000412_min_max_vals = [
        [int(e.split("\t")[1]), int(e.split("\t")[2])]
        for e in expected_sorted_exclusions
        if e.startswith("ENST00000000412")
    ]

    with open(paths.exclusions_shuffled, "r") as f:
        shuffled_content = f.readlines()

    for line in shuffled_content:
        line = line.strip()
        chrom, start, stop, *other = line.split("\t")
        if chrom == "ENST00000000412":
            start = int(start)
            stop = int(stop)
            for input_start, input_stop in ENST00000000412_min_max_vals:
                overlap = (start <= input_stop) and (input_start <= stop)
                assert not overlap


@pytest.mark.dependency(depends=["filter_trans_files"])
def test__randomize_with_target_file(test_files):
    paths, transcripts, auxiliaries = test_files

    prep_coords.filter_transcript_files(paths, transcripts)

    with open(paths.filtered_transcript, "r") as f:
        expected_trans_content = f.readlines()
        expected_trans_content.append(tab_line("ENST00000000233", 543))

    with open(paths.filtered_protein_transcript, "r") as f:
        expected_trans_prot_content = f.readlines()
        expected_trans_prot_content.append(tab_line("ENST00000000233", 180))

    prep_coords._randomize_with_target_file(paths, transcripts, seed=1234)

    check_expected_content(expected_trans_content, paths.filtered_transcript)
    check_expected_content(
        expected_trans_prot_content, paths.filtered_protein_transcript
    )

    # TODO: Complete


def test__non_randomized(test_files):
    paths, transcripts, auxiliaries = test_files
    # expected_content = [
    #     tab_line("ENST00000000233", 115, 124),
    #     tab_line("ENST00000000233", 164, 177),
    #     tab_line("ENST00000000412", 113, 124),
    #     tab_line("ENST00000001008", 27, 36),
    #     tab_line("ENST00000001008", 189, 198),
    # ]
    prep_coords._non_randomized(paths)

    # TODO: See docstring in _non_randomized
    # check_expected_content(expected_content, paths.exclusions_shuffled)


def test__exclude_positively_selected_genes_disabled(test_files):
    paths, transcripts, auxiliaries = test_files

    # Dummy data written to shuffled exclusions file:
    # When positively selected genes are disabled, should just copy this file
    paths.exclusions_shuffled.write_text(tab_line("chr1", 123, 456))

    prep_coords._exclude_positively_selected_genes_disabled(paths)

    with open(paths.exclusions_shuffled, "r") as f:
        expected_content = f.readlines()

    check_expected_content(expected_content, paths.epitopes)


def test__exclude_positively_selected_genes(test_files):
    paths, transcripts, auxiliaries = test_files

    # Write dummy data for "shuffle file"
    paths.exclusions_shuffled.write_text(
        tab_line("ENST00000000233", 164, 177)
        + tab_line("ENST00000001008", 113, 124)
    )

    # In the exclusions aux file, we have ENST00000001008
    expected_content = [tab_line("ENST00000000233", 164, 177)]
    prep_coords._exclude_positively_selected_genes(paths, auxiliaries)

    check_expected_content(expected_content, paths.epitopes)


def test_get_protein_complement(minimal_epitopes):
    # TODO: undo this fixture implementation... more explicit using
    #       standard test_files

    paths = minimal_epitopes

    # The expected complement file
    expected_tmp_content = [
        tab_line("ENST00000000233", 0, 100),
        tab_line("ENST00000000233", 150, 400),
        tab_line("ENST00000000412", 0, 277),
        tab_line("ENST00000001008", 0, 459),
    ]

    # NOTE: the length in the filtered protein must be >= the upperbound of
    # the largest stop position in the epitopes file
    # (for the associated chroms), otherwise
    # bedtools will raise a warning which is undetected via the subprocess

    # NOTE: In the absence of start, stop in the epitopes file,
    # chromosomes appearing in the filtered file will automatically be
    # provided the complementary interval [0, length]

    prep_coords.get_protein_complement(paths)
    check_expected_content(expected_tmp_content, paths.intra_epitopes_prot_tmp)

    # The final step of the calculation essentially finds word regex matches
    # with the initial input. Therefore, in this case, the expected content are
    # the lines starting with ENST00000000233
    check_expected_content(expected_tmp_content[:2], paths.intra_epitopes_prot)


def test__prep_ssb192(test_files):
    paths, *others = test_files

    start_stop_233 = ("ENST00000000233", 400, 500)
    start_stop_234 = ("ENST00000000233", 100, 150)
    start_stop_235 = ("ENST00000000235", 0, 10)

    paths.epitopes.write_text(
        tab_line(*start_stop_233)
        + tab_line(*start_stop_234)
        + tab_line(*start_stop_235)
    )

    def _map(start_stop: tuple):
        chrom, start, stop = start_stop

        c1 = start * 3 - 6
        if c1 < 0:
            c1 += 3
        c2 = stop * 3 + 3

        return chrom, c1, c2, *start_stop

    expected_content = [
        tab_line(*_map(start_stop_233)),
        tab_line(*_map(start_stop_234)),
        tab_line(*_map(start_stop_235)),
    ]

    prep_coords._prep_ssb192(paths)

    check_expected_content(expected_content, paths.epitopes_cds)


def test__prep_not_ssb192(test_files):
    paths, *other = test_files

    start_stop_233 = ("ENST00000000233", 400, 500)
    start_stop_234 = ("ENST00000000234", 100, 150)
    start_stop_235 = ("ENST00000000235", 0, 10)

    paths.epitopes.write_text(
        tab_line(*start_stop_233)
        + tab_line(*start_stop_234)
        + tab_line(*start_stop_235)
    )

    prep_coords._prep_not_ssb192(paths)

    def _map(start_stop: tuple):
        chrom, start, stop = start_stop

        c1 = start * 3 - 3
        c2 = stop * 3

        return chrom, c1, c2, *start_stop

    expected_content = [
        tab_line(*_map(start_stop_233)),
        tab_line(*_map(start_stop_234)),
        tab_line(*_map(start_stop_235)),
    ]

    check_expected_content(expected_content, paths.epitopes_cds)


def test_transform_coordinates(test_files):
    paths, *others = test_files

    start_stop_233 = ("ENST00000000233", 400, 500)
    start_stop_234 = ("ENST00000000234", 100, 150)
    start_stop_235 = ("ENST00000000235", 10, 20)

    # Build input epitopes file
    paths.epitopes.write_text(
        tab_line(*start_stop_233)
        + tab_line(*start_stop_234)
        + tab_line(*start_stop_235)
    )

    def _map(start_stop: tuple):
        chrom, start, stop = start_stop

        c1 = start * 3 - 3
        c2 = stop * 3

        return chrom, c1, c2, *start_stop

    # Build input epitopes cds coordinate file
    paths.epitopes_cds.write_text(
        tab_line(*_map(start_stop_233))
        + tab_line(*_map(start_stop_234))
        + tab_line(*_map(start_stop_235))
    )
    # Which should be ...
    # ENST00000000233 1197    1500    ENST00000000233 400     500
    # ENST00000000234 297     450     ENST00000000234 100     150
    # ENST00000000235 27      60      ENST00000000235 10      20

    # Some compatible transcript content
    paths.filtered_transcript.write_text(
        tab_line("ENST00000000233", 1500)
        + tab_line("ENST00000000234", 500)
        + tab_line("ENST00000000235", 80)
    )

    expected_content_tmp = [
        tab_line("ENST00000000233", 0, 1197),
        tab_line("ENST00000000234", 0, 297),
        tab_line("ENST00000000234", 450, 500),
        tab_line("ENST00000000235", 0, 27),
        tab_line("ENST00000000235", 60, 80),
    ]
    prep_coords.transform_coordinates(paths)
    check_expected_content(expected_content_tmp, paths.intra_epitopes_tmp)

    check_expected_content(expected_content_tmp, paths.intra_epitopes)
