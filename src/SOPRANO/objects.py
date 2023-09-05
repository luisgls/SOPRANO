import pathlib
from argparse import Namespace
from dataclasses import dataclass


def _data_dir():
    return pathlib.Path(__file__).parent.joinpath("data")


@dataclass(frozen=True)
class TranscriptPaths:
    transcript_length: pathlib.Path
    protein_transcript_length: pathlib.Path


EnsemblTranscripts = TranscriptPaths(
    transcript_length=_data_dir().joinpath("ensemble_transcript.length"),
    protein_transcript_length=_data_dir().joinpath(
        "ensemble_transcript_protein.length"
    ),
)


@dataclass(frozen=True)
class _GenomicPaths:
    sizes: pathlib.Path
    fasta: pathlib.Path


GRCh37 = _GenomicPaths(
    sizes=_data_dir().joinpath("chrom_GRCh37.sizes"),
    fasta=_data_dir().joinpath(""),  # TODO: implement
)

GRCh38 = _GenomicPaths(
    sizes=_data_dir().joinpath("chrom_GRCh38.sizes"),  # TODO: implement
    fasta=_data_dir().joinpath(""),  # TODO: implement
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
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        target_regions: pathlib.Path | None = None,
    ):
        self.analysis_name = analysis_name
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
    "use_ssb192",
    "use_random",
    "exclude_drivers",
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
        transcripts: TranscriptPaths,
    ):
        super().__init__(analysis_name, bed_path, cache_dir, target_regions)

        self.transcripts = transcripts
        self.use_ssb192 = use_ssb192
        self.use_random = use_random
        self.exclude_drivers = exclude_drivers

    @classmethod
    def from_namespace(cls, namespace: Namespace):
        for k in namespace.__dict__.keys():
            assert k in _NAMESPACE_KEYS, k

        transcripts = TranscriptPaths(
            namespace.transcript, namespace.protein_transcript
        )

        return cls(
            namespace.analysis_name,
            namespace.input_path,
            namespace.bed_path,
            namespace.cache_dir,
            namespace.target_regions,
            namespace.use_ssb192,
            namespace.use_random,
            namespace.exclude_drivers,
            transcripts,
        )
