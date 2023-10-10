import pathlib
from argparse import Namespace
from dataclasses import dataclass
from typing import Set

from SOPRANO.utils.path_utils import Directories, genome_pars_to_paths
from SOPRANO.utils.url_utils import (
    build_ensembl_urls,
    check_ensembl_file_url,
    compute_chrom_sizes,
    compute_fasta_index,
    decompress,
    download_from_url,
    filename_from_url,
    find_earliest_release,
    find_latest_release,
)


@dataclass(frozen=True)
class TranscriptPaths:
    transcript_length: pathlib.Path
    protein_transcript_length: pathlib.Path
    transcript_fasta: pathlib.Path

    @classmethod
    def defaults(cls):
        return cls(
            transcript_length=Directories.data("ensemble_transcript.length"),
            protein_transcript_length=Directories.data(
                "ensemble_transcript_protein.length"
            ),
            transcript_fasta=Directories.data("ensemble_transcriptID.fasta"),
        )


@dataclass(frozen=True)
class GenomePaths:
    sizes: pathlib.Path
    fasta: pathlib.Path

    @classmethod
    def GRCh37(cls, release=110):
        fasta, sizes = genome_pars_to_paths("GRCh37", release)
        return cls(sizes=sizes, fasta=fasta)

    @classmethod
    def GRCh38(cls, release=110):
        fasta, sizes = genome_pars_to_paths("GRCh38", release)
        return cls(sizes=sizes, fasta=fasta)


@dataclass(frozen=True)
class AuxiliaryPaths:
    genes_to_exclude: pathlib.Path
    intron_length: pathlib.Path

    @classmethod
    def defaults(cls):
        return cls(
            genes_to_exclude=Directories.data("genes2exclude.txt"),
            intron_length=Directories.data("transcript_intron_length.bed"),
        )


class AnalysisPaths:
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        random_regions: pathlib.Path | None = None,
    ):
        self.analysis_name = analysis_name
        self.input_path = input_path
        self.bed_path = bed_path
        self.random_regions_path = random_regions
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
        file_name = f"{self.analysis_name}"

        if len(extensions) > 0:
            file_name += "." + ".".join([*extensions])

        return self.cache_dir.joinpath(file_name)


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
    "keep_drivers",
    "seed",
    "species",
    "assembly",
    "release",
)


class Parameters(AnalysisPaths):
    def __init__(
        self,
        analysis_name: str,
        input_path: pathlib.Path,
        bed_path: pathlib.Path,
        cache_dir: pathlib.Path,
        random_regions: pathlib.Path | None,
        use_ssb192: bool,
        use_random: bool,
        exclude_drivers: bool,
        seed: int,
        transcripts: TranscriptPaths,
        genomes: GenomePaths,
    ):
        super().__init__(
            analysis_name, input_path, bed_path, cache_dir, random_regions
        )

        self.transcripts = transcripts
        self.genomes = genomes
        self.use_ssb192 = use_ssb192
        self.use_random_regions = random_regions is not None
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

        species = namespace.species
        assembly = namespace.assembly
        release = namespace.release

        if assembly == "GRCh37":
            assert species == "homo_sapiens", species
            genomes = (
                EnsemblData.homo_sapiens_GRCh37().get_genome_reference_paths(
                    release
                )
            )
        else:
            genomes = EnsemblData(
                species=species, assembly=assembly
            ).get_genome_reference_paths(release)

        return cls(
            analysis_name=namespace.analysis_name,
            input_path=namespace.input_path,
            bed_path=namespace.bed_path,
            cache_dir=namespace.cache_dir,
            random_regions=namespace.random_regions,
            use_ssb192=namespace.use_ssb192,
            use_random=namespace.use_random,
            exclude_drivers=not namespace.keep_drivers,
            seed=namespace.seed,
            transcripts=transcripts,
            genomes=genomes,
        )


class SOPRANOError(Exception):
    pass


class _GatherReferences:
    # Urls
    toplevel_url: str
    primary_assembly_url: str

    # Params
    species: str
    assembly: str

    # Status
    toplevel_gz_done: Set[int] = set()
    toplevel_fa_done: Set[int] = set()
    toplevel_fai_done: Set[int] = set()
    primary_assembly_gz_done: Set[int] = set()
    primary_assembly_fa_done: Set[int] = set()
    primary_assembly_fai_done: Set[int] = set()
    sizes_done: Set[int] = set()

    def _dest_directory(self, release: int):
        return Directories.data(self.species) / f"{release}_{self.assembly}"

    def _dest_fa_gz(self, release: int, _toplevel: bool):
        return self._dest_directory(release) / filename_from_url(
            self.toplevel_url if _toplevel else self.primary_assembly_url
        ).format(RELEASE=release)

    def _dest_fa(self, release: int, _toplevel: bool):
        return self._dest_fa_gz(release, _toplevel).with_suffix("")

    def _dest_fai(self, release: int, _toplevel: bool):
        return self._dest_fa_gz(release, _toplevel).with_suffix(".fai")

    def dest_chrom(self, release: int, _toplevel: bool):
        return self._dest_fa(release, _toplevel).with_suffix(".chrom")

    def toplevel_fa_gz_path(self, release: int):
        return self._dest_fa_gz(release, _toplevel=True)

    def toplevel_fa_path(self, release: int):
        return self._dest_fa(release, _toplevel=True)

    def toplevel_fai_path(self, release: int):
        return self._dest_fai(release, _toplevel=True)

    def toplevel_chrom_path(self, release: int):
        return self.toplevel_fa_path(release).with_suffix(".chrom")

    def primary_assembly_fa_gz_path(self, release: int):
        return self._dest_fa_gz(release, _toplevel=False)

    def primary_assembly_fa_path(self, release: int):
        return self._dest_fa(release, _toplevel=False)

    def primary_assembly_fai_path(self, release: int):
        return self._dest_fai(release, _toplevel=False)

    def _download(self, release: int, _toplevel):
        if _toplevel:
            source_url = self.toplevel_url
            dest_path = self.toplevel_fa_gz_path(release)
            decompressed_path = self.toplevel_fa_path(release)
        else:
            source_url = self.primary_assembly_url
            dest_path = self.primary_assembly_fa_gz_path(release)
            decompressed_path = self.primary_assembly_fa_path(release)

        if not (decompressed_path.exists() or dest_path.exists()):
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            check_ensembl_file_url(source_url, release)
            download_from_url(
                source_url.format(RELEASE=release),
                target_path=dest_path,
            )

    def _check_release_ok(self, release):
        min_release = find_earliest_release(self.toplevel_url)
        max_release = find_latest_release(self.toplevel_url)

        if not (min_release <= release <= max_release):
            raise ValueError(release)

    def download_toplevel(self, release):
        if release not in self.toplevel_gz_done:
            self._check_release_ok(release)

            if not self.toplevel_fa_gz_path(release).exists():
                self._download(release, _toplevel=True)

            self.toplevel_gz_done.add(release)

    def download_primary_assembly(self, release):
        if release not in self.primary_assembly_gz_done:
            self._check_release_ok(release)

            if not self.primary_assembly_fa_gz_path(release).exists():
                self._download(release, _toplevel=False)

            self.primary_assembly_gz_done.add(release)

    def decompress_toplevel(self, release):
        if release not in self.toplevel_fa_done:
            if not self.toplevel_fa_path(release).exists():
                decompress(self.toplevel_fa_gz_path(release))

            self.toplevel_fa_done.add(release)

    def decompress_primary_assembly(self, release):
        if release not in self.primary_assembly_fa_done:
            if not self.primary_assembly_fa_path(release).exists():
                decompress(self.primary_assembly_fa_gz_path(release))

            self.primary_assembly_fa_done.add(release)

    def compute_chrom_sizes(self, release):
        if release not in self.sizes_done:
            if not self.toplevel_chrom_path(release).exists():
                compute_chrom_sizes(self.toplevel_fai_path(release))

            self.sizes_done.add(release)

    def compute_fasta_index_toplevel(self, release):
        if release not in self.toplevel_fai_done:
            if not self.toplevel_fai_path(release).exists():
                compute_fasta_index(self.toplevel_fa_path(release))

            self.toplevel_fai_done.add(release)

    def compute_fasta_index_primary_assembly(self, release):
        if release not in self.primary_assembly_fai_done:
            if not self.primary_assembly_fai_path(release).exists():
                compute_fasta_index(self.primary_assembly_fa_path(release))

            self.primary_assembly_fai_done.add(release)

    def compute_all_toplevel(self, release):
        self.download_toplevel(release)
        self.decompress_toplevel(release)
        self.compute_fasta_index_toplevel(release)
        self.compute_chrom_sizes(release)

    def compute_all_primary_assembly(self, release):
        self.download_primary_assembly(release)
        self.decompress_primary_assembly(release)
        self.compute_fasta_index_primary_assembly(release)

    def get_genome_reference_paths(self, release):
        return GenomePaths(
            sizes=self.toplevel_chrom_path(release),
            fasta=self.toplevel_fa_path(release),
        )


class EnsemblData(_GatherReferences):
    def __init__(self, species: str, assembly: str, _init_urls=True):
        self.species = species
        self.assembly = assembly

        if _init_urls:
            url_dict = build_ensembl_urls(species, assembly)
            self.toplevel_url = url_dict["toplevel"]
            self.primary_assembly_url = url_dict["primary_assembly"]

    @classmethod
    def homo_sapiens_GRCh38(cls):
        return cls("homo_sapiens", "GRCh38")

    @classmethod
    def homo_sapiens_GRCh37(cls):
        # GRCh37 is actually has a deviant url structure, so manually set here
        toplevel_url = (
            "https://ftp.ensembl.org/pub/grch37/release-{RELEASE}/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.toplevel.fa.gz"
        )

        primary_assembly_url = (
            "https://ftp.ensembl.org/pub/grch37/release-{RELEASE}/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
        )

        species = "homo_sapiens"
        assembly = "GRCh37"

        obj = cls(species, assembly, _init_urls=False)
        obj.toplevel_url = toplevel_url
        obj.primary_assembly_url = primary_assembly_url
        return obj
