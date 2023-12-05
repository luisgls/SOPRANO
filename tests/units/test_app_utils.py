import os
from pathlib import Path

import pytest

from SOPRANO.core.objects import GenomePaths
from SOPRANO.utils.app_utils import (
    DownloaderUIOptions,
    DownloaderUIProcessing,
    ImmunopeptidomesUIOptions,
    ImmunopeptidomeUIProcessing,
    LinkVEPUIProcessing,
    PipelineUIOptions,
    PipelineUIProcessing,
    _lines_ok,
    _select_from_dict,
    process_text_and_file_inputs,
)
from SOPRANO.utils.path_utils import _SOPRANO_DEFAULT_CACHE, Directories


class TemporaryFiles:
    def __init__(self, pass_ext, fail_ext, target_dir, mkdir, files):
        self._pass_ext = pass_ext
        self._fail_ext = fail_ext
        self._target_dir = target_dir
        self._mkdir = mkdir
        self._files = files

    @classmethod
    def tmp_in_dir(cls, target_dir: Path, ext: str):
        files = target_dir / f"pytest.{ext}", target_dir / "pytest.fail"
        return cls(
            pass_ext=files[0],
            fail_ext=files[1],
            target_dir=target_dir,
            mkdir=False,
            files=files,
        )

    @classmethod
    def tmp_files(cls, target_dir: Path, *names: str, mkdir=False):
        files = [target_dir / n for n in names]
        return cls(
            pass_ext=None,
            fail_ext=None,
            target_dir=target_dir,
            mkdir=mkdir,
            files=files,
        )

    def __enter__(self):
        if self._mkdir:
            self._target_dir.mkdir(exist_ok=False)

        for f in self._files:
            f.touch(exist_ok=False)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for f in self._files:
            f.unlink(missing_ok=False)

        if self._mkdir:
            self._target_dir.rmdir()


def _check_options_generator(method, directory: Path, extension: str):
    with TemporaryFiles.tmp_in_dir(directory, extension) as c:
        options = method()
        assert c._pass_ext.name in options
        assert options[c._pass_ext.name] == c._pass_ext
        assert c._fail_ext.name not in options
        with pytest.raises(KeyError):
            assert options[c._fail_ext.name] == c._fail_ext


def test__select_from_dict():
    k, v = "foo", "bar"
    assert _select_from_dict(k, {k: v}) == v


def test_pipeline_options_genome_reference(mock_genome_dir):
    assembly, release, mock_dir = mock_genome_dir

    with TemporaryFiles.tmp_files(
        mock_dir,
        "something.dna.toplevel.fa",
        "something.dna.toplevel.chrom",
        mkdir=True,
    ):
        options = PipelineUIOptions.genome_reference()
        assert "{} - Ensembl release {}".format(assembly, release) in options


def test_pipeline_options_annotated_mutations():
    _check_options_generator(
        PipelineUIOptions.annotated_mutations,
        Directories.app_annotated_inputs(),
        "anno",
    )


def test_pipeline_options_immunopeptidome():
    _check_options_generator(
        PipelineUIOptions.immunopeptidome,
        Directories.app_immunopeptidomes(),
        "bed",
    )


def test_pipeline_options_substitution_method():
    for v in (192, 7):
        assert v in PipelineUIOptions.substitution_method().values()


def test_pipeline_options_coordinates():
    _check_options_generator(
        PipelineUIOptions.coordinates,
        Directories.app_coordinate_files(),
        "bed",
    )


def test_pipeline_processing_genome_reference(mock_genome_dir):
    assembly, release, mock_dir = mock_genome_dir

    selection_string = "{} - Ensembl release {}".format(assembly, release)

    condition, output = PipelineUIProcessing.genome_reference(selection_string)

    assert isinstance(output, GenomePaths)
    assert output.fasta.name.endswith(".fa")
    assert output.sizes.name.endswith(".chrom")


def test_pipeline_processing_annotated_mutations():
    with TemporaryFiles.tmp_in_dir(
        Directories.app_annotated_inputs(), "anno"
    ) as c:
        condition, processed_option = PipelineUIProcessing.annotated_mutations(
            c._pass_ext.name
        )

        assert processed_option == c._pass_ext

        with pytest.raises(KeyError):
            PipelineUIProcessing.annotated_mutations(c._fail_ext.name)


def test_pipeline_processing_immunopeptidome():
    with TemporaryFiles.tmp_in_dir(
        Directories.app_immunopeptidomes(), "bed"
    ) as c:
        (
            condition,
            processed_immunopeptidome,
        ) = PipelineUIProcessing.immunopeptidome(c._pass_ext.name)

        assert processed_immunopeptidome == c._pass_ext

        with pytest.raises(KeyError):
            PipelineUIProcessing.immunopeptidome(c._fail_ext.name)


def test_pipeline_processing_substitution_method():
    condition, processed_192 = PipelineUIProcessing.substitution_method(
        "SSB192"
    )
    assert processed_192 == 192
    condition, processed_7 = PipelineUIProcessing.substitution_method("SSB7")
    assert processed_7 == 7

    with pytest.raises(KeyError):
        PipelineUIProcessing.substitution_method("SSB000")


def test_pipeline_processing_job_name():
    cache_key = "SOPRANO_CACHE"
    if cache_key in os.environ:
        _soprano_cache = os.environ[cache_key]
        del os.environ[cache_key]
    else:
        _soprano_cache = _SOPRANO_DEFAULT_CACHE.as_posix()
    try:
        mock_dir = Path("/just/for/pytest")
        os.environ[cache_key] = mock_dir.as_posix()
        job_name = "soprano"
        condition, processed_name = PipelineUIProcessing.job_name(job_name)
        assert processed_name == mock_dir / job_name
    finally:
        os.environ[cache_key] = _soprano_cache


def test_vep_processing_cache_location():
    test_loc = Path.cwd()
    assert LinkVEPUIProcessing.cache_location(test_loc.as_posix()) == test_loc


def test_downloader_options_type():
    options = DownloaderUIOptions.type()
    assert "toplevel" in options
    assert "primary_assembly" in options


def test_downloader_processing_species():
    assert DownloaderUIProcessing.species("Big Dog") == "big_dog"


def test_downloader_processing_assembly():
    assert DownloaderUIProcessing.assembly("cat") == "cat"


def test_downloader_processing_release():
    assert DownloaderUIProcessing.release("123") == 123


def test_downloader_processing_type():
    with pytest.raises(ValueError):
        DownloaderUIProcessing.type("cat")

    for permitted in ("primary_assembly", "toplevel"):
        assert DownloaderUIProcessing.type(permitted) == permitted


def test_immunopeptidome_options_hla_alleles():
    options = ImmunopeptidomesUIOptions.hla_alleles()

    assert isinstance(options, list)
    assert len(options) == 354  # NOTE: Based on hard coded file... 13 Oct 2023


def test_immunopeptidome_options_transcript_ids():
    options = ImmunopeptidomesUIOptions.transcript_ids()

    assert isinstance(options, list)
    assert (
        len(options) == 9754
    )  # NOTE: Based on hard coded file... 13 Oct 2023


def test_immunopeptidome_options_subset_method():
    expected = {"None", "Retention", "Exclusion"}

    assert set(ImmunopeptidomesUIOptions.subset_method()) == expected


def test_immunopeptidome_processing_hla_alleles():
    assert ImmunopeptidomeUIProcessing.hla_alleles([1, 2, 3]) == [1, 2, 3]


def test_immunopeptidome_processing_transcript_ids():
    assert ImmunopeptidomeUIProcessing.transcript_ids([1, 2, 3]) == [1, 2, 3]


def test_immunopeptidome_processing_subset_method():
    assert ImmunopeptidomeUIProcessing.subset_method([], "foo") == ([], [])
    assert ImmunopeptidomeUIProcessing.subset_method([1, 2, 3], "None") == (
        [],
        [],
    )

    assert ImmunopeptidomeUIProcessing.subset_method([1, 2], "Retention") == (
        [1, 2],
        [],
    )
    assert ImmunopeptidomeUIProcessing.subset_method([1, 2], "Exclusion") == (
        [],
        [1, 2],
    )
    with pytest.raises(ValueError):
        ImmunopeptidomeUIProcessing.subset_method([1, 2], "bar")


def test_immunopeptidome_processing_name():
    assert ImmunopeptidomeUIProcessing.name("x") == "x.bed"
    assert ImmunopeptidomeUIProcessing.name("x.bed") == "x.bed"


def test__lines_ok():
    assert _lines_ok(["a", "b", "c"], 0, 10)
    assert not _lines_ok(["a", "b", "c"], 0, 2)
    with pytest.raises(ValueError):
        _lines_ok(["a", "b", "c"], 4, 3)


def test_process_text_and_file_inputs():
    assert process_text_and_file_inputs("cats\nand\ndogs") == (
        True,
        ["cats", "and", "dogs"],
    )
