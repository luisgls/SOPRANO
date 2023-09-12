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


GRCh37 = GenomePaths(
    sizes=_data_dir().joinpath("chrom_GRCh37.sizes"),
    fasta=_data_dir().joinpath("Homo_sapiens.GRCh37.dna.toplevel.fa"),
)

GRCh38 = GenomePaths(
    sizes=_data_dir().joinpath("chrom_GRCh38.sizes"),  # TODO: implement
    fasta=_data_dir().joinpath("Homo_sapiens.GRCh38.dna.toplevel.fa"),
)


@dataclass(frozen=True)
class AuxiliaryPaths:
    genes_to_exclude: pathlib.Path


AuxiliaryFiles = AuxiliaryPaths(
    genes_to_exclude=_data_dir().joinpath("genes2exclude.txt")
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
        self.variants_silent = self._cached_path("variants", "silent")
        self.variants_nonsilent = self._cached_path("variants", "nonsilent")
        self.variants_missense = self._cached_path(
            "variants", "variants_missense"
        )

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
)


class Parameters(AnalysisPaths):
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        target_regions: pathlib.Path,
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
        assert set(namespace.__dict__.keys()) == set(
            _NAMESPACE_KEYS
        ), namespace.__dict__.keys()

        transcripts = TranscriptPaths(
            namespace.transcript,
            namespace.protein_transcript,
            namespace.transcript_ids,
        )

        if namespace.genome_ref == "grch37":
            genomes = GRCh37
        elif namespace.genome_ref == "grch38":
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
