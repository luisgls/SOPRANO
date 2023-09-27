import pathlib
from argparse import Namespace
from dataclasses import dataclass


def _data_dir():
    return pathlib.Path(__file__).parent.joinpath("data")


@dataclass(frozen=True)
class TranscriptPaths:
    transcript_length: pathlib.Path
    protein_transcript_length: pathlib.Path
    transcript_fasta: pathlib.Path


EnsemblTranscripts = TranscriptPaths(
    transcript_length=_data_dir().joinpath("ensemble_transcript.length"),
    protein_transcript_length=_data_dir().joinpath(
        "ensemble_transcript_protein.length"
    ),
    transcript_fasta=_data_dir().joinpath("ensemble_transcriptID.fasta"),
)


@dataclass(frozen=True)
class GenomePaths:
    sizes: pathlib.Path
    fasta: pathlib.Path


def _genome_pars_to_paths(ref, release):
    data_dir = (
        _data_dir().joinpath("homo_sapiens").joinpath(f"{release}_{ref}")
    )

    genome_path = data_dir.joinpath(f"Homo_sapiens.{ref}.dna.toplevel.fa")
    chroms_path = data_dir.joinpath(f"Homo_sapiens.{ref}.dna.toplevel.chrom")
    return genome_path, chroms_path


GRCh37 = GenomePaths(
    sizes=_genome_pars_to_paths("GRCh37", 110)[1],
    fasta=_genome_pars_to_paths("GRCh37", 110)[0],
)

GRCh38 = GenomePaths(
    sizes=_genome_pars_to_paths("GRCh38", 110)[1],
    fasta=_genome_pars_to_paths("GRCh38", 110)[0],
)


@dataclass(frozen=True)
class AuxiliaryPaths:
    genes_to_exclude: pathlib.Path
    intron_length: pathlib.Path


AuxiliaryFiles = AuxiliaryPaths(
    genes_to_exclude=_data_dir().joinpath("genes2exclude.txt"),
    intron_length=_data_dir().joinpath("transcript_intron_length.bed"),
)


def cache_path_builder(tmpdir: pathlib.Path, name: str, *extensions: str):
    file_name = f"{name}"

    if len(extensions) > 0:
        file_name += "." + ".".join([*extensions])

    return tmpdir.joinpath(file_name)


class AnalysisPaths:
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        target_regions: pathlib.Path | None = None,
    ):
        self.analysis_name = analysis_name
        self.input_path = input_path
        self.bed_path = bed_path
        self.target_regions_path = target_regions
        self.cache_dir = cache_dir

        # Transcripts
        self.filtered_protein_transcript = self._cached_path(
            "protein_length_filt", "txt"
        )
        self.filtered_transcript = self._cached_path(
            "transcript_length_filt", "txt"
        )

        # Exclusion regions for randomization
        self.exclusions = self._cached_path("exclusion", "ori")
        self.exclusions_sorted = self._cached_path("exclusion", "bed")
        self.exclusions_shuffled = self._cached_path("epitopes", "ori2")

        # Epitope files
        self.epitopes = self._cached_path("epitopes", "bed")
        self.epitopes_cds = self._cached_path("epitopes_cds", "bed")
        self.epitopes_cds_fasta = self._cached_path("epitopes_cds", "fasta")
        self.epitopes_trans_regs = self._cached_path(
            "epitopes_cds", "transcript_regions"
        )
        self.epitopes_trans_regs_sum = self._cached_path(
            "epitopes_cds", "transcript_regions", "sum"
        )

        # Complement (intra) epitope files
        self.intra_epitopes = self._cached_path("intra_epitopes", "bed")
        self.intra_epitopes_tmp = self._cached_path(
            "intra_epitopes", "bed", "tmp"
        )
        self.intra_epitopes_prot = self._cached_path(
            "intra_epitopes_prot", "bed"
        )
        self.intra_epitopes_prot_tmp = self._cached_path(
            "intra_epitopes_prot", "bed", "tmp"
        )
        self.intra_epitopes_cds = self._cached_path("intra_epitopes_cds")
        self.intra_epitopes_cds_fasta = self._cached_path(
            "intra_epitopes_cds", "fasta"
        )
        self.intra_epitopes_trans_regs = self._cached_path(
            "intra_epitopes_cds", "transcript_regions"
        )
        self.intra_epitopes_trans_regs_sum = self._cached_path(
            "intra_epitopes_cds", "transcript_regions", "sum"
        )

        # Analysis files
        self.sim_fixed = self._cached_path("sim_fixed")
        self.col_corrected = self._cached_path("col_corrected")
        self.contextualised = self._cached_path("contextualised")
        self.flagged = self._cached_path("flagged")
        self.triplet_counts = self._cached_path("triplets", "counts")
        self.final_epitope_corrections = self._cached_path(
            "corrected_matrix", "epitopes"
        )
        self.final_intra_epitope_corrections = self._cached_path(
            "corrected_matrix", "intra_epitopes"
        )
        self.epitope_nans = self._cached_path("epitopes", "nans")
        self.intra_epitope_nans = self._cached_path("intra_epitopes", "nans")

        # variant counts
        self.variants_silent = self._cached_path("variants", "silent", "bed")
        self.variants_nonsilent = self._cached_path(
            "variants", "nonsilent", "bed"
        )
        self.variants_missense = self._cached_path(
            "variants", "missense", "bed"
        )
        self.variants_intronic = self._cached_path(
            "variants", "intronic", "bed"
        )
        self.raw_silent_count = self._cached_path("raw", "silent", "count")
        self.raw_nonsilent_count = self._cached_path(
            "raw", "nonsilent", "count"
        )
        self.raw_missense_count = self._cached_path("raw", "missense", "count")
        self.in_silent_count = self._cached_path("in", "silent", "count")
        self.in_nonsilent_count = self._cached_path("in", "nonsilent", "count")
        self.in_missense_count = self._cached_path("in", "missense", "count")
        self.out_silent_count = self._cached_path("out", "silent", "count")
        self.out_nonsilent_count = self._cached_path(
            "out", "nonsilent", "count"
        )
        self.out_missense_count = self._cached_path("out", "missense", "count")

        self.data_epitopes = self._cached_path("data", "epitopes")

        self.intron_rate = self._cached_path("intron", "rate")

        self.results_path = self._cached_path("results.tsv")

    def _cached_path(self, *extensions):
        return cache_path_builder(
            self.cache_dir, self.analysis_name, *extensions
        )


_NAMESPACE_KEYS = (
    "analysis_name",
    "input_path",
    "bed_path",
    "cache_dir",
    "target_regions",
    "transcript",
    "protein_transcript",
    "transcript_ids",
    "use_ssb192",
    "use_random",
    "exclude_drivers",
    "seed",
    "genome_ref",
    "release",
)


class Parameters(AnalysisPaths):
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        target_regions: pathlib.Path | None,
        use_ssb192: bool,
        use_random: bool,
        exclude_drivers: bool,
        seed: int,
        transcripts: TranscriptPaths,
        genomes: GenomePaths,
    ):
        super().__init__(
            analysis_name, input_path, bed_path, cache_dir, target_regions
        )

        self.transcripts = transcripts
        self.genomes = genomes
        self.use_ssb192 = use_ssb192
        self.use_target_regions = target_regions is not None
        self.use_random = use_random
        self.exclude_drivers = exclude_drivers
        self.seed = None if seed < 0 else seed

    @classmethod
    def from_namespace(cls, namespace: Namespace):
        for k in namespace.__dict__.keys():
            assert k in _NAMESPACE_KEYS, k

        for k in _NAMESPACE_KEYS:
            assert k in namespace.__dict__.keys(), k

        transcripts = TranscriptPaths(
            namespace.transcript,
            namespace.protein_transcript,
            namespace.transcript_ids,
        )

        if namespace.genome_ref == "GRCh37":
            genomes = GRCh37
        elif namespace.genome_ref == "GRCh38":
            genomes = GRCh38
        else:
            raise KeyError(f"Unrecognized reference: {namespace.genome_ref}")

        return cls(
            namespace.analysis_name,
            namespace.input_path,
            namespace.bed_path,
            namespace.cache_dir,
            namespace.target_regions,
            namespace.use_ssb192,
            namespace.use_random,
            namespace.exclude_drivers,
            namespace.seed,
            transcripts,
            genomes,
        )
