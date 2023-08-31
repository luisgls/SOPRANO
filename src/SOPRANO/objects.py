import pathlib
from dataclasses import dataclass


def _data_dir():
    return pathlib.Path(__file__).parent.joinpath("data")


@dataclass(frozen=True)
class TranscriptPaths:
    transcript_length: pathlib.Path
    protein_transcript_length: pathlib.Path


@dataclass(frozen=True)
class _GenomicPaths:
    sizes: pathlib.Path
    fasta: pathlib.Path


EnsemblTranscripts = TranscriptPaths(
    transcript_length=_data_dir().joinpath("ensemble_transcript.length"),
    protein_transcript_length=_data_dir().joinpath(
        "ensemble_transcript_protein.length"
    ),
)

GRCh37 = _GenomicPaths(
    sizes=_data_dir().joinpath("chrom_GRCh37.sizes"),
    fasta=_data_dir().joinpath(""),  # TODO: implement
)

GRCh38 = _GenomicPaths(
    sizes=_data_dir().joinpath("chrom_GRCh38.sizes"),  # TODO: implement
    fasta=_data_dir().joinpath(""),  # TODO: implement
)


@dataclass(frozen=True)
class _AuxiliaryPaths:
    genes_to_exclude: pathlib.Path


AuxiliaryPaths = _AuxiliaryPaths(
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
        tmpdir: pathlib.Path,
        target_regions_path: pathlib.Path | None = None,
    ):
        self.name = analysis_name
        self.bed_path = bed_path
        self.target_regions_path = target_regions_path
        self.tmpdir = tmpdir

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
        self.intra_epitopes_prot = self._cached_path(
            "intra_epitopes_prot", "bed"
        )

    def _cached_path(self, *extensions):
        return cache_path_builder(self.tmpdir, self.name, *extensions)
