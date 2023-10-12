from contextlib import contextmanager, redirect_stdout
from io import StringIO

import streamlit as st

from SOPRANO.core.objects import EnsemblData
from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe


@contextmanager
def st_capture(output_func):
    """Taken from https://discuss.streamlit.io/t/
        cannot-print-the-terminal-output-in-streamlit/6602/2

    :param output_func: an instance of st.empty()

    """
    with StringIO() as stdout, redirect_stdout(stdout):
        old_write = stdout.write

        def new_write(string):
            ret = old_write(string)
            output_func(stdout.getvalue())
            return ret

        stdout.write = new_write
        yield


def _select_from_dict(selection: str, selection_dict: dict):
    selection_value = selection_dict[selection]
    st.text(f"Selected: {selection_value}")
    return selection_value


class _PipelineUI:
    @staticmethod
    def genome_reference(*args, **kwargs):
        pass

    @staticmethod
    def annotated_mutations(*args, **kwargs):
        pass

    @staticmethod
    def immunopeptidome(*args, **kwargs):
        pass

    @staticmethod
    def substitution_method(*args, **kwargs):
        pass

    @staticmethod
    def coordinates(*args, **kwargs):
        pass


class PipelineUIOptions(_PipelineUI):
    @staticmethod
    def genome_reference():
        homo_sapiens_dir = Directories.genomes_homo_sapiens()

        genome_dirs = [
            item for item in homo_sapiens_dir.glob("*") if item.is_dir()
        ]

        # Remove bad options (i.e. no toplevel fa and chrom files)
        for item in genome_dirs[::-1]:
            toplevel_path = item.glob("*dna*toplevel*.fa")
            chrom_path = item.glob("*dna*toplevel*.chrom")

            if len(list(toplevel_path)) == len(list(chrom_path)) == 1:
                pass
            else:
                genome_dirs.remove(item)

        genome_ids = [
            "{} - Ensembl release {}".format(*x.name.split("_")[::-1])
            for x in genome_dirs
        ]

        options_dict = {
            name: dir_path for name, dir_path in zip(genome_ids, genome_dirs)
        }

        return options_dict

    @staticmethod
    def annotated_mutations():
        options_dict = {}
        for directory in (
            Directories.examples(),
            Directories.app_annotated_inputs(),
        ):
            for x in directory.glob("*.anno*"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def immunopeptidome():
        options_dict = {}
        for directory in (
            Directories.immunopeptidomes_humans(),
            Directories.app_immunopeptidomes(),
        ):
            for x in directory.glob("*.bed"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def substitution_method():
        return {"SSB192": 192, "SSB7": 7}

    @staticmethod
    def coordinates():
        options_dict = {None: None}
        for x in Directories.app_coordinate_files().glob("*.bed"):
            options_dict[x.name] = x
        return options_dict


class PipelineUIProcessing(_PipelineUI):
    @staticmethod
    def genome_reference(genome_selection: str):
        assembly, release = genome_selection.split(" - Ensembl release ")
        data = EnsemblData(species="homo_sapiens", assembly=assembly)
        fasta_path = data.toplevel_fa_path(int(release))
        chrom_path = data.toplevel_chrom_path(int(release))
        st.text(f"Selected: {fasta_path}, {chrom_path}")
        return data.get_genome_reference_paths(int(release))

    @staticmethod
    def annotated_mutations(annotation_selection: str):
        options_dict = PipelineUIOptions.annotated_mutations()
        return _select_from_dict(annotation_selection, options_dict)

    @staticmethod
    def immunopeptidome(immunopeptidome_selection: str):
        options_dict = PipelineUIOptions.immunopeptidome()
        return _select_from_dict(immunopeptidome_selection, options_dict)

    @staticmethod
    def substitution_method(subs_selection: str):
        options_dict = PipelineUIOptions.substitution_method()
        return _select_from_dict(subs_selection, options_dict)

    @staticmethod
    def coordinates(coordinates_selection: str):
        options_dict = PipelineUIOptions.coordinates()
        return _select_from_dict(coordinates_selection, options_dict)

    @staticmethod
    def job_name(job_name: str):
        cache_dir = Directories.cache(job_name)
        st.text(f"Selected: {cache_dir}")
        return cache_dir


class _LinkVEPUI:
    @staticmethod
    def cache_location(*args, **kwargs):
        pass


class LinkVEPUIOptions(_LinkVEPUI):
    pass


class LinkVEPUIProcessing(_LinkVEPUI):
    @classmethod
    def cache_location(cls):
        pass


class _DownloaderUI:
    @staticmethod
    def species(*args, **kwargs):
        pass

    @staticmethod
    def assembly(*args, **kwargs):
        pass

    @staticmethod
    def release(*args, **kwargs):
        pass

    @staticmethod
    def type(*args, **kwargs):
        pass


class DownloaderUIOptions(_DownloaderUI):
    @staticmethod
    def type():
        return {"toplevel": "toplevel", "primary_assembly": "primary_assembly"}


class DownloaderUIProcessing(_DownloaderUI):
    @staticmethod
    def species(species: str):
        return species.lower().replace(" ", "_")

    @staticmethod
    def type(type_selection: str):
        return type_selection


class _AnnotatorUI:
    pass


class AnnotatorUIOptions(_AnnotatorUI):
    pass


class AnnotatorUIProcessing(_AnnotatorUI):
    pass


class _ImmunopeptidomeUI:
    @staticmethod
    def hla_alleles(*args, **kwargs):
        pass

    @staticmethod
    def transcript_ids(*args, **kwargs):
        pass


class ImmunopeptidomeUIProcessing(_ImmunopeptidomeUI):
    pass


class ImmunopeptidomesUIOptions(_ImmunopeptidomeUI):
    @staticmethod
    def hla_alleles():
        hla_types_path = Directories.examples("TCGA_hlaTypesAll.tsv")
        options = pipe(
            ["cut", "-f3", hla_types_path.as_posix()],
            ["tr", ",", "\n"],
            ["sort", "-u"],
        ).split("\n")
        return options

    @staticmethod
    def transcript_ids():
        hla_binders_path = Directories.data(
            "allhlaBinders_exprmean1.IEDBpeps.bed.unique_ids"
        )

        with open(hla_binders_path, "r") as f:
            transcript_options = f.read()

        return transcript_options.split("\n")
