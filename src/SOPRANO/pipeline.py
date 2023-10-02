from typing import List

from SOPRANO.core.analysis import (
    _build_flag_file,
    _check_triplet_counts,
    _col_correct,
    _compute_theoretical_subs_192,
    _context_correction,
    _correct_from_total_sites,
    _fix_simulated,
    _initial_triplet_counts,
    _sum_possible_across_region,
)
from SOPRANO.core.dnds import _compute_coverage, _intersect_introns
from SOPRANO.core.intersect import (
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
from SOPRANO.core.objects import AuxiliaryPaths, Parameters
from SOPRANO.core.obtain_fasta_regions import (
    _get_non_target_regions,
    _get_target_fasta_regions,
    _get_trans_regs,
)
from SOPRANO.core.prepare_coordinates import (
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
from SOPRANO.utils.path_utils import _check_paths
from SOPRANO.utils.print_utils import task_output

AuxiliaryFiles = AuxiliaryPaths.defaults()


class _PipelineComponent:
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

    def check_ready(self, params: Parameters):
        _check_paths(
            params.filtered_transcript,
            params.filtered_protein_transcript,
        )


class NonRandom(_Randomize):
    msg = "No randomization selected; selecting unique items from bed file"

    def _apply(self, params: Parameters):
        _non_randomized(params)


class RandomizeWithoutRegions(_Randomize):
    msg = "Performing randomization without supplement bed file definitions"

    def _apply(self, params: Parameters):
        _define_excluded_regions_for_randomization(params)
        _sort_excluded_regions_for_randomization(params, seed=params.seed)


class RandomizeWithRegions(_Randomize):
    msg = "Performing randomization using supplement bed file definitions"

    def _apply(self, params: Parameters):
        _randomize_with_target_file(
            params, params.transcripts, seed=params.seed
        )


class _GeneExclusions(_PipelineComponent):
    def check_ready(self, params: Parameters):
        _check_paths(params.exclusions_shuffled)


class GeneExclusions(_GeneExclusions):
    msg = "Excluding positively selected genes"

    def _apply(self, params: Parameters):
        _exclude_positively_selected_genes(params, AuxiliaryFiles)


class GeneExclusionsDisabled(_GeneExclusions):
    msg = "Retaining positively selected genes"

    def _apply(self, params: Parameters):
        _exclude_positively_selected_genes_disabled(params)


class BuildProteinComplement(_PipelineComponent):
    msg = "Building protein complement"

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes, params.filtered_protein_transcript)

    def _apply(self, params: Parameters):
        _get_protein_complement(params)


class _SSB192Selection(_PipelineComponent):
    """Intermediate class for ssb192 mutrate selection"""

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes)


class UseSSB192(_SSB192Selection):
    msg = "Preparing CDS coordinates using SSB 192 substitution"

    def _apply(self, params: Parameters):
        _prep_ssb192(params)


class NotSSB192(_SSB192Selection):
    msg = "Preparing CDS coordinates using SSB 7 substitution"

    def _apply(self, params: Parameters):
        _prep_not_ssb192(params)


class BuildIntraEpitopesCDS(_PipelineComponent):
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
    msg = "Compiling list of transcript:regions to estimate number of sites"

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes_cds_fasta, params.intra_epitopes_cds)

    def _apply(self, params: Parameters):
        _get_trans_regs(params.epitopes_cds_fasta, params.epitopes_trans_regs)
        _get_trans_regs(
            params.intra_epitopes_cds_fasta, params.intra_epitopes_trans_regs
        )


class ComputeSSB192TheoreticalSubs(_PipelineComponent):
    msg = "Computing all theoretical substitutions for SSB192"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.epitopes_cds_fasta,
            params.intra_epitopes_cds_fasta,
            params.epitopes_trans_regs,
            params.intra_epitopes_trans_regs,
        )

    def _apply(self, params: Parameters):
        _compute_theoretical_subs_192(
            params.epitopes_cds_fasta, params.epitopes_trans_regs
        )
        _compute_theoretical_subs_192(
            params.intra_epitopes_cds_fasta, params.intra_epitopes_trans_regs
        )


class SumPossibleAcrossRegions(_PipelineComponent):
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
    msg = "Processing VEP annotated file to estimated 192 rate parameters"
    # TODO: Fix msg

    def check_ready(self, params: Parameters):
        # TODO
        pass

    def _apply(self, params: Parameters):
        _fix_simulated(params)


class ColumnCorrect(_PipelineComponent):
    msg = "Applying column corrections"

    def check_ready(self, params: Parameters):
        _check_paths(params.sim_fixed)

    def _apply(self, params: Parameters):
        _col_correct(params)


class ContextCorrection(_PipelineComponent):
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
    msg = "Flagging calculations"

    def check_ready(self, params: Parameters):
        _check_paths(params.col_corrected, params.contextualised)

    def _apply(self, params: Parameters):
        _build_flag_file(params)


class TripletCounts(_PipelineComponent):
    msg = "Counting triplets"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.col_corrected, params.contextualised, params.flagged
        )

    def _apply(self, params: Parameters):
        _initial_triplet_counts(params)


class SiteCorrections(_PipelineComponent):
    msg = "Performing site corrections"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.epitopes_trans_regs_sum,
            params.intra_epitopes_trans_regs_sum,
            params.triplet_counts,
        )
        _check_triplet_counts(params)

    def _apply(self, params: Parameters):
        _correct_from_total_sites(params)


class IntersectByFrequency(_PipelineComponent):
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
    msg = "Computing silent mutation counts"

    def _apply(self, params: Parameters):
        _get_silent_variant_counts(params)


class GetNonSilentCounts(_PipelineComponent):
    msg = "Computing non-silent mutation counts"

    def _apply(self, params: Parameters):
        _get_nonsilent_variant_counts(params)


class GetMissenseCounts(_PipelineComponent):
    msg = "Computing missense mutation counts"

    def _apply(self, params: Parameters):
        _get_missense_variant_counts(params)


class GetIntronicCounts(_PipelineComponent):
    msg = "Computing intronic mutation counts"

    def _apply(self, params: Parameters):
        _get_intronic_variant_counts(params)


class OnOffCounts(_PipelineComponent):
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
    msg = "Computing intronic rate"

    def check_ready(self, params: Parameters):
        _check_paths(params.variants_intronic)

    def _apply(self, params: Parameters):
        _intersect_introns(params, AuxiliaryFiles)


class ComputeStatistics(_PipelineComponent):
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
    jobs: List[_PipelineComponent] = [FilterTranscripts()]

    if params.use_random_regions:
        jobs.append(RandomizeWithRegions())
    elif params.use_random:
        jobs.append(RandomizeWithoutRegions())
    else:
        jobs.append(NonRandom())

    if params.exclude_drivers:
        jobs.append(GeneExclusions())
    else:
        jobs.append(GeneExclusionsDisabled())

    jobs.append(BuildProteinComplement())

    if params.use_ssb192:
        jobs.append(UseSSB192())
    else:
        jobs.append(NotSSB192())

    jobs.append(BuildIntraEpitopesCDS())
    jobs.append(ObtainFastaRegions())
    jobs.append(GetTranscriptRegionsForSites())

    if params.use_ssb192:
        jobs.append(ComputeSSB192TheoreticalSubs())
        jobs.append(SumPossibleAcrossRegions())
        jobs.append(FixSimulated())
        jobs.append(ColumnCorrect())
        jobs.append(ContextCorrection())
        jobs.append(FlagComputations())
        jobs.append(TripletCounts())
        jobs.append(SiteCorrections())
    else:
        raise KeyError("SSB7 requires implementation")

    jobs.append(IntersectByFrequency())
    jobs.append(GetSilentCounts())
    jobs.append(GetNonSilentCounts())
    jobs.append(GetMissenseCounts())
    jobs.append(GetIntronicCounts())
    jobs.append(OnOffCounts())
    jobs.append(BuildEpitopesDataFile())
    jobs.append(CheckTargetMutations())
    jobs.append(ComputeIntronRate())
    jobs.append(ComputeStatistics())

    for job in jobs:
        job.apply(params)
