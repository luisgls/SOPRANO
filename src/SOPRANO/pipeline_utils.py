import pathlib
from datetime import datetime
from typing import List

from SOPRANO.analysis import (
    _build_flag_file,
    _check_triplet_counts,
    _col_correct,
    _compute_theoretical_subs,
    _context_correction,
    _correct_from_total_sites,
    _fix_simulated,
    _initial_triplet_counts,
    _sum_possible_across_region,
)
from SOPRANO.dnds import _compute_coverage, _intersect_introns
from SOPRANO.intersect import (
    _check_target_mutations,
    _count_intersected_mutations,
    _count_mutations,
    _get_intronic_variant_counts,
    _get_missense_variant_counts,
    _get_nonsilent_variant_counts,
    _get_silent_variant_counts,
    _intersect_by_frequency,
    _update_epitopes_data_file,
)
from SOPRANO.misc_utils import MissingDataError
from SOPRANO.objects import AuxiliaryFiles, Parameters
from SOPRANO.obtain_fasta_regions import (
    _get_non_target_regions,
    _get_target_fasta_regions,
    _get_trans_regs,
)
from SOPRANO.prepare_coordinates import (
    _define_excluded_regions_for_randomization,
    _exclude_positively_selected_genes,
    _exclude_positively_selected_genes_disabled,
    _get_protein_complement,
    _non_randomized,
    _prep_not_ssb192,
    _prep_ssb192,
    _randomize_with_target_file,
    _sort_excluded_regions_for_randomization,
    filter_transcript_files,
    transform_coordinates,
)


def time_output():
    now = datetime.now()
    return now.strftime("%d/%m/%Y %H:%M:%S")


def task_output(msg):
    print(f"[{time_output()}] {msg}")


def _check_paths(*dependent_paths: pathlib.Path):
    for path in dependent_paths:
        if not path.exists():
            raise MissingDataError(path)


class _PipelineComponent:
    """
    Components of the pipeline are designed to follow the pattern:

    Component.apply(params)

    where apply should include the call to check_read()
    to permit execution.

    Pipeline components should override these methods.
    """

    @staticmethod
    def apply(params: Parameters):
        pass

    @staticmethod
    def check_ready(params: Parameters):
        pass


class _PipelineComponent2:
    msg: str | None = None  # Will be printed to stdout with date/time stamp

    def apply(self, params: Parameters):
        if self.msg is not None:
            task_output(self.msg)
        self.check_ready(params)
        self._apply(params)

    def _apply(self, params: Parameters):
        pass

    def check_ready(self, params: Parameters):
        pass


class FilterTranscripts(_PipelineComponent):
    """
    Filter transcript files with respect to input bed file
    """

    @staticmethod
    def apply(params: Parameters):
        FilterTranscripts.check_ready(params)
        filter_transcript_files(params, params.transcripts)

    @staticmethod
    def check_ready(params: Parameters):
        for path in (
            params.bed_path,
            params.transcripts.transcript_length,
            params.transcripts.protein_transcript_length,
        ):
            if not path.exists():
                raise MissingDataError(path)


class FilterTranscripts2(_PipelineComponent2):
    msg = "Filtering transcripts"

    def _apply(self, params: Parameters):
        filter_transcript_files(params, params.transcripts)

    def check_ready(self, params: Parameters):
        _check_paths(
            params.bed_path,
            params.transcripts.transcript_length,
            params.transcripts.protein_transcript_length,
        )


class _Randomize(_PipelineComponent):
    """Intermediate class for randomization procedures"""

    @staticmethod
    def check_ready(params: Parameters):
        for path in (
            params.filtered_transcript,
            params.filtered_protein_transcript,
        ):
            if not path.exists():
                raise MissingDataError(
                    f"Filtered transcript not found: {path}"
                )


class _Randomize2(_PipelineComponent2):
    """Intermediate class for randomization procedures"""

    def check_ready(self, params: Parameters):
        _check_paths(
            params.filtered_transcript,
            params.filtered_protein_transcript,
        )


class NonRandom(_Randomize):
    """No randomization implemented"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _non_randomized(params)


class NonRandom2(_Randomize2):
    msg = "No randomization selected; selecting unique items from bed file"

    def _apply(self, params: Parameters):
        _non_randomized(params)


class RandomizeWithoutRegions(_Randomize):
    """Randomizes without user input file"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _define_excluded_regions_for_randomization(params)
        _sort_excluded_regions_for_randomization(params, seed=params.seed)


class RandomizeWithoutRegions2(_Randomize2):
    msg = "Performing randomization without supplement bed file definitions"

    def _apply(self, params: Parameters):
        _define_excluded_regions_for_randomization(params)
        _sort_excluded_regions_for_randomization(params, seed=params.seed)


class RandomizeWithRegions(_Randomize):
    """Randomizes with user input file"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _randomize_with_target_file(
            params, params.transcripts, seed=params.seed
        )


class RandomizeWithRegions2(_Randomize2):
    msg = "Performing randomization using supplement bed file definitions"

    def _apply(self, params: Parameters):
        _randomize_with_target_file(
            params, params.transcripts, seed=params.seed
        )


class _GeneExclusions(_PipelineComponent):
    """Intermediate class for gene exclusions"""

    @staticmethod
    def check_ready(params: Parameters):
        if not params.exclusions_shuffled.exists():
            raise MissingDataError(f"{params.exclusions_shuffled}")


class _GeneExclusions2(_PipelineComponent2):
    def check_ready(self, params: Parameters):
        _check_paths(params.exclusions_shuffled)


class GeneExclusions(_GeneExclusions):
    """Applies gene exclusions"""

    @staticmethod
    def apply(params: Parameters):
        _GeneExclusions.check_ready(params)
        _exclude_positively_selected_genes(params, AuxiliaryFiles)


class GeneExclusions2(_GeneExclusions2):
    msg = "Excluding positively selected genes"

    def _apply(self, params: Parameters):
        _exclude_positively_selected_genes(params, AuxiliaryFiles)


class GeneExclusionsDisabled(_GeneExclusions):
    """No gene exclusions"""

    @staticmethod
    def apply(params: Parameters):
        _GeneExclusions.check_ready(params)
        _exclude_positively_selected_genes_disabled(params)


class GeneExclusionsDisabled2(_GeneExclusions2):
    msg = "Retaining positively selected genes"

    def _apply(self, params: Parameters):
        _exclude_positively_selected_genes_disabled(params)


class BuildProteinComplement(_PipelineComponent):
    """Build protein complement file"""

    @staticmethod
    def check_ready(params: Parameters):
        for path in (params.epitopes, params.filtered_protein_transcript):
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        BuildProteinComplement.check_ready(params)
        _get_protein_complement(params)


class BuildProteinComplement2(_PipelineComponent2):
    msg = "Building protein complement"

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes, params.filtered_protein_transcript)

    def _apply(self, params: Parameters):
        _get_protein_complement(params)


class _SSB192Selection(_PipelineComponent):
    """Intermediate class for ssb192 mutrate selection"""

    @staticmethod
    def check_ready(params: Parameters):
        if not params.epitopes.exists():
            raise MissingDataError(params.epitopes.as_posix())


class _SSB192Selection2(_PipelineComponent2):
    """Intermediate class for ssb192 mutrate selection"""

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes)


class UseSSB192(_SSB192Selection):
    """Applies ssb192 selection in CDS coordinate prep"""

    @staticmethod
    def apply(params: Parameters):
        _SSB192Selection.check_ready(params)
        _prep_ssb192(params)


class UseSSB1922(_SSB192Selection2):
    msg = "Preparing CDS coordinates using SSB 192 substitution"

    def _apply(self, params: Parameters):
        _prep_ssb192(params)


class NotSSB192(_SSB192Selection):
    """Does not apply ssb192 selection in CDS coordiante prep"""

    @staticmethod
    def apply(params: Parameters):
        _SSB192Selection.check_ready(params)
        _prep_not_ssb192(params)


class NotSSB1922(_SSB192Selection2):
    msg = "Preparing CDS coordinates using SSB 7 substitution"

    def _apply(self, params: Parameters):
        _prep_not_ssb192(params)


class BuildIntraEpitopesCDS(_PipelineComponent):
    """Builds the complement to the epitope in cds coords"""

    @staticmethod
    def check_ready(params: Parameters):
        for path in (
            params.epitopes,
            params.epitopes_cds,
            params.filtered_transcript,
        ):
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        BuildIntraEpitopesCDS.check_ready(params)
        transform_coordinates(params)


class BuildIntraEpitopesCDS2(_PipelineComponent2):
    msg = "Building intra epitope file in CDS coordinates"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.epitopes,
            params.epitopes_cds,
            params.filtered_transcript,
        )

    def _apply(self, params: Parameters):
        transform_coordinates(params)


class ObtainFastaRegions(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.transcripts.transcript_fasta,
            params.epitopes_cds,
            params.intra_epitopes_cds,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ObtainFastaRegions.check_ready(params)
        _get_target_fasta_regions(params, params.transcripts)
        _get_non_target_regions(params, params.transcripts)


class ObtainFastaRegions2(_PipelineComponent2):
    msg = "Obtaining fasta regions"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.transcripts.transcript_fasta,
            params.epitopes_cds,
            params.intra_epitopes_cds,
        )

    def _apply(self, params: Parameters):
        _get_target_fasta_regions(params, params.transcripts)
        _get_non_target_regions(params, params.transcripts)


class GetTranscriptRegionsForSites(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (params.epitopes_cds_fasta, params.intra_epitopes_cds)
        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        GetTranscriptRegionsForSites.check_ready(params)
        _get_trans_regs(params.epitopes_cds_fasta, params.epitopes_trans_regs)
        _get_trans_regs(
            params.intra_epitopes_cds_fasta, params.intra_epitopes_trans_regs
        )


class GetTranscriptRegionsForSites2(_PipelineComponent2):
    msg = "Compiling list of transcript:regions to estimate number of sites"

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes_cds_fasta, params.intra_epitopes_cds)

    def _apply(self, params: Parameters):
        _get_trans_regs(params.epitopes_cds_fasta, params.epitopes_trans_regs)
        _get_trans_regs(
            params.intra_epitopes_cds_fasta, params.intra_epitopes_trans_regs
        )


class ComputeSSB192TheoreticalSubs(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.epitopes_cds_fasta,
            params.intra_epitopes_cds_fasta,
            params.epitopes_trans_regs,
            params.intra_epitopes_trans_regs,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ComputeSSB192TheoreticalSubs.check_ready(params)
        _compute_theoretical_subs(
            params.epitopes_cds_fasta, params.epitopes_trans_regs
        )
        _compute_theoretical_subs(
            params.intra_epitopes_cds_fasta, params.intra_epitopes_trans_regs
        )


class ComputeSSB192TheoreticalSubs2(_PipelineComponent2):
    msg = "Computing all theoretical substitutions for SSB192"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.epitopes_cds_fasta,
            params.intra_epitopes_cds_fasta,
            params.epitopes_trans_regs,
            params.intra_epitopes_trans_regs,
        )

    def _apply(self, params: Parameters):
        _compute_theoretical_subs(
            params.epitopes_cds_fasta, params.epitopes_trans_regs
        )
        _compute_theoretical_subs(
            params.intra_epitopes_cds_fasta, params.intra_epitopes_trans_regs
        )


class SumPossibleAcrossRegions(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (params.epitopes_trans_regs, params.intra_epitopes_trans_regs)
        for p in paths:
            if not p.exists():
                raise MissingDataError(p)

    @staticmethod
    def apply(params: Parameters):
        SumPossibleAcrossRegions.check_ready(params)
        _sum_possible_across_region(
            params.epitopes_trans_regs, params.epitopes_trans_regs_sum
        )
        _sum_possible_across_region(
            params.intra_epitopes_trans_regs,
            params.intra_epitopes_trans_regs_sum,
        )


class SumPossibleAcrossRegions2(_PipelineComponent2):
    msg = "Computing sum over possible sites in on and off target regions"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.epitopes_trans_regs, params.intra_epitopes_trans_regs
        )

    def _apply(self, params: Parameters):
        _sum_possible_across_region(
            params.epitopes_trans_regs, params.epitopes_trans_regs_sum
        )
        _sum_possible_across_region(
            params.intra_epitopes_trans_regs,
            params.intra_epitopes_trans_regs_sum,
        )


class FixSimulated(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        FixSimulated.check_ready(params)
        _fix_simulated(params)


class FixSimulated2(_PipelineComponent2):
    msg = "Processing VEP annotated file to estimated 192 rate parameters"
    # TODO: Fix msg

    def check_ready(self, params: Parameters):
        # TODO
        pass

    def _apply(self, params: Parameters):
        _fix_simulated(params)


class ColumnCorrect(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        if not params.sim_fixed.exists():
            raise MissingDataError(params.sim_fixed)

    @staticmethod
    def apply(params: Parameters):
        ColumnCorrect.check_ready(params)
        _col_correct(params)


class ColumnCorrect2(_PipelineComponent2):
    msg = "Applying column corrections"

    def check_ready(self, params: Parameters):
        _check_paths(params.sim_fixed)

    def _apply(self, params: Parameters):
        _col_correct(params)


class ContextCorrection(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.col_corrected,
            params.genomes.fasta,
            params.genomes.sizes,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ContextCorrection.check_ready(params)
        _context_correction(params, params.genomes)


class ContextCorrection2(_PipelineComponent2):
    msg = "Applying context corrections"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.col_corrected,
            params.genomes.fasta,
            params.genomes.sizes,
        )

    def _apply(self, params: Parameters):
        _context_correction(params, params.genomes)


class FlagComputations(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (params.col_corrected, params.contextualised)
        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        FlagComputations.check_ready(params)
        _build_flag_file(params)


class FlagComputations2(_PipelineComponent2):
    msg = "Flagging calculations"

    def check_ready(self, params: Parameters):
        _check_paths(params.col_corrected, params.contextualised)

    def _apply(self, params: Parameters):
        _build_flag_file(params)


class TripletCounts(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (params.col_corrected, params.contextualised, params.flagged)
        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        TripletCounts.check_ready(params)
        _initial_triplet_counts(params)


class TripletCounts2(_PipelineComponent2):
    msg = "Counting triplets"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.col_corrected, params.contextualised, params.flagged
        )

    def _apply(self, params: Parameters):
        _initial_triplet_counts(params)


class SiteCorrections(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.epitopes_trans_regs_sum,
            params.intra_epitopes_trans_regs_sum,
            params.triplet_counts,
        )

        for p in paths:
            if not p.exists():
                raise MissingDataError(p)

        _check_triplet_counts(params)

    @staticmethod
    def apply(params: Parameters):
        SiteCorrections.check_ready(params)
        _correct_from_total_sites(params)


class SiteCorrections2(_PipelineComponent2):
    msg = "Performing site corrections"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.epitopes_trans_regs_sum,
            params.intra_epitopes_trans_regs_sum,
            params.triplet_counts,
        )

    def _apply(self, params: Parameters):
        _correct_from_total_sites(params)


class IntersectByFrequency(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.final_epitope_corrections,
            params.final_intra_epitope_corrections,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        IntersectByFrequency.check_ready(params)
        _intersect_by_frequency(
            params.final_epitope_corrections, params.epitope_nans
        )
        _intersect_by_frequency(
            params.final_intra_epitope_corrections, params.intra_epitope_nans
        )


class IntersectByFrequency2(_PipelineComponent2):
    msg = "Intersecting by frequencies"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.final_epitope_corrections,
            params.final_intra_epitope_corrections,
        )

    def _apply(self, params: Parameters):
        _intersect_by_frequency(
            params.final_epitope_corrections, params.epitope_nans
        )
        _intersect_by_frequency(
            params.final_intra_epitope_corrections, params.intra_epitope_nans
        )


class GetSilentCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_silent_variant_counts(params)


class GetSilentCounts2(_PipelineComponent2):
    msg = "Computing silent mutation counts"

    def _apply(self, params: Parameters):
        _get_silent_variant_counts(params)


class GetNonSilentCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_nonsilent_variant_counts(params)


class GetNonSilentCounts2(_PipelineComponent2):
    msg = "Computing non-silent mutation counts"

    def _apply(self, params: Parameters):
        _get_nonsilent_variant_counts(params)


class GetMissenseCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_missense_variant_counts(params)


class GetMissenseCounts2(_PipelineComponent2):
    msg = "Computing missense mutation counts"

    def _apply(self, params: Parameters):
        _get_missense_variant_counts(params)


class GetIntronicCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_intronic_variant_counts(params)


class GetIntronicCounts2(_PipelineComponent2):
    msg = "Computing intronic mutation counts"

    def _apply(self, params: Parameters):
        _get_intronic_variant_counts(params)


class OnOffCounts(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.variants_silent,
            params.variants_nonsilent,
            params.variants_missense,
            params.variants_intronic,
            params.epitopes,
            params.intra_epitopes_prot,
        )

        for p in paths:
            if not p.exists():
                raise MissingDataError(p)

    @staticmethod
    def apply(params: Parameters):
        raw_silent = _count_mutations(
            params.variants_silent, params.raw_silent_count
        )
        raw_nonsilent = _count_mutations(
            params.variants_nonsilent, params.raw_nonsilent_count
        )
        raw_missense = _count_mutations(
            params.variants_missense, params.raw_missense_count
        )
        in_silent = _count_intersected_mutations(
            params.variants_silent, params.epitopes, params.in_silent_count
        )
        in_nonsilent = _count_intersected_mutations(
            params.variants_nonsilent,
            params.epitopes,
            params.in_nonsilent_count,
        )
        in_missense = _count_intersected_mutations(
            params.variants_missense, params.epitopes, params.in_missense_count
        )
        out_silent = _count_intersected_mutations(
            params.variants_silent,
            params.intra_epitopes_prot,
            params.out_silent_count,
        )
        out_nonsilent = _count_intersected_mutations(
            params.variants_nonsilent,
            params.intra_epitopes_prot,
            params.out_nonsilent_count,
        )
        out_missense = _count_intersected_mutations(
            params.variants_missense,
            params.intra_epitopes_prot,
            params.out_missense_count,
        )

        def _print(region, silent, nonsilent, missense):
            print(f"{region}:")
            print("{0:.<30}".format("Silent") + silent)
            print("{0:.<30}".format("Non-silent") + nonsilent)
            print("{0:.<30}".format("Missense") + missense)

        _print("Global region", raw_silent, raw_nonsilent, raw_missense)
        _print("(ON) Target region", in_silent, in_nonsilent, in_missense)
        _print("(OFF) Target region", out_silent, out_nonsilent, out_missense)


class OnOffCounts2(_PipelineComponent2):
    msg = "Computing summary of global/on/off region counts"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.variants_silent,
            params.variants_nonsilent,
            params.variants_missense,
            params.variants_intronic,
            params.epitopes,
            params.intra_epitopes_prot,
        )

    def _apply(self, params: Parameters):
        raw_silent = _count_mutations(
            params.variants_silent, params.raw_silent_count
        )
        raw_nonsilent = _count_mutations(
            params.variants_nonsilent, params.raw_nonsilent_count
        )
        raw_missense = _count_mutations(
            params.variants_missense, params.raw_missense_count
        )
        in_silent = _count_intersected_mutations(
            params.variants_silent, params.epitopes, params.in_silent_count
        )
        in_nonsilent = _count_intersected_mutations(
            params.variants_nonsilent,
            params.epitopes,
            params.in_nonsilent_count,
        )
        in_missense = _count_intersected_mutations(
            params.variants_missense, params.epitopes, params.in_missense_count
        )
        out_silent = _count_intersected_mutations(
            params.variants_silent,
            params.intra_epitopes_prot,
            params.out_silent_count,
        )
        out_nonsilent = _count_intersected_mutations(
            params.variants_nonsilent,
            params.intra_epitopes_prot,
            params.out_nonsilent_count,
        )
        out_missense = _count_intersected_mutations(
            params.variants_missense,
            params.intra_epitopes_prot,
            params.out_missense_count,
        )

        def _print(region, silent, nonsilent, missense):
            print(f"{region}:")
            print("{0:.<30}".format("Silent") + silent)
            print("{0:.<30}".format("Non-silent") + nonsilent)
            print("{0:.<30}".format("Missense") + missense)

        _print("Global region", raw_silent, raw_nonsilent, raw_missense)
        _print("(ON) Target region", in_silent, in_nonsilent, in_missense)
        _print("(OFF) Target region", out_silent, out_nonsilent, out_missense)


class BuildEpitopesDataFile(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.variants_silent,
            params.variants_nonsilent,
            params.in_silent_count,
            params.in_nonsilent_count,
            params.out_silent_count,
            params.out_nonsilent_count,
        )
        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        BuildEpitopesDataFile.check_ready(params)

        variants = (params.variants_silent, params.variants_nonsilent)
        in_counts = (params.in_silent_count, params.in_nonsilent_count)
        out_counts = (params.out_silent_count, params.out_nonsilent_count)
        labs = ("synonymous", "missense")

        for variant_count, in_out_count, use_epi, lab in zip(
            variants * 2,
            in_counts + out_counts,
            (True, True, False, False),
            labs * 2,
        ):
            _update_epitopes_data_file(
                variant_count,
                in_out_count,
                params,
                _use_epitope=use_epi,
                _label=lab,
            )


class BuildEpitopesDataFile2(_PipelineComponent2):
    msg = "Building epitope combined data file"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.variants_silent,
            params.variants_nonsilent,
            params.in_silent_count,
            params.in_nonsilent_count,
            params.out_silent_count,
            params.out_nonsilent_count,
        )

    def _apply(self, params: Parameters):
        variants = (params.variants_silent, params.variants_nonsilent)
        in_counts = (params.in_silent_count, params.in_nonsilent_count)
        out_counts = (params.out_silent_count, params.out_nonsilent_count)
        labs = ("synonymous", "missense")

        for variant_count, in_out_count, use_epi, lab in zip(
            variants * 2,
            in_counts + out_counts,
            (True, True, False, False),
            labs * 2,
        ):
            _update_epitopes_data_file(
                variant_count,
                in_out_count,
                params,
                _use_epitope=use_epi,
                _label=lab,
            )


class CheckTargetMutations(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.in_silent_count,
            params.in_nonsilent_count,
            params.in_missense_count,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        CheckTargetMutations.check_ready(params)
        _check_target_mutations(params)


class CheckTargetMutations2(_PipelineComponent2):
    msg = "Checking target mutations"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.in_silent_count,
            params.in_nonsilent_count,
            params.in_missense_count,
        )

    def _apply(self, params: Parameters):
        _check_target_mutations(params)


class ComputeIntronRate(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        if not params.variants_intronic.exists():
            raise MissingDataError(params.variants_intronic)

    @staticmethod
    def apply(params: Parameters):
        _intersect_introns(params, AuxiliaryFiles)


class ComputeIntronRate2(_PipelineComponent2):
    msg = "Computing intronic rate"

    def check_ready(self, params: Parameters):
        _check_paths(params.variants_intronic)

    def _apply(self, params: Parameters):
        _intersect_introns(params, AuxiliaryFiles)


class ComputeStatistics(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.data_epitopes,
            params.epitope_nans,
            params.intra_epitope_nans,
            params.intron_rate,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ComputeStatistics.check_ready(params)
        _compute_coverage(params)


class ComputeStatistics2(_PipelineComponent2):
    msg = "Computing dN/dS statistical summary"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.data_epitopes,
            params.epitope_nans,
            params.intra_epitope_nans,
            params.intron_rate,
        )

    def _apply(self, params: Parameters):
        _compute_coverage(params)


def run_pipeline(params: Parameters):
    jobs: List[_PipelineComponent2] = [FilterTranscripts2()]

    if params.use_target_regions:
        jobs.append(RandomizeWithRegions2())
    elif params.use_random:
        jobs.append(RandomizeWithoutRegions2())
    else:
        jobs.append(NonRandom2())

    if params.exclude_drivers:
        jobs.append(GeneExclusions2())
    else:
        jobs.append(GeneExclusionsDisabled2())

    jobs.append(BuildProteinComplement2())

    if params.use_ssb192:
        jobs.append(UseSSB1922())
    else:
        jobs.append(NotSSB1922())

    jobs.append(BuildIntraEpitopesCDS2())
    jobs.append(ObtainFastaRegions2())
    jobs.append(GetTranscriptRegionsForSites2())

    if params.use_ssb192:
        jobs.append(ComputeSSB192TheoreticalSubs2())
        jobs.append(SumPossibleAcrossRegions2())
        jobs.append(FixSimulated2())
        jobs.append(ColumnCorrect2())
        jobs.append(ContextCorrection2())
        jobs.append(FlagComputations2())
        jobs.append(TripletCounts2())
        jobs.append(SiteCorrections2())
    else:
        raise KeyError("SSB7 requires implementation")

    jobs.append(IntersectByFrequency2())
    jobs.append(GetSilentCounts2())
    jobs.append(GetNonSilentCounts2())
    jobs.append(GetMissenseCounts2())
    jobs.append(GetIntronicCounts2())
    jobs.append(OnOffCounts2())
    jobs.append(BuildEpitopesDataFile2())
    jobs.append(CheckTargetMutations2())
    jobs.append(ComputeIntronRate2())
    jobs.append(ComputeStatistics2())

    for job in jobs:
        job.apply(params)
