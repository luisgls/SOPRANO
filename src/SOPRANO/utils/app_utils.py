from contextlib import contextmanager, redirect_stdout
from io import StringIO

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


class AppOptions:
    @staticmethod
    def get_human_genome_options():
        homo_sapiens_dir = Directories.genomes_homo_sapiens()

        genome_dirs = [
            item for item in homo_sapiens_dir.glob("*") if item.is_dir()
        ]

        # Remove unviable options (i.e. no toplevel fa and chrom files)
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
    def get_annotated_input_options():
        options_dict = {}
        for directory in (
            Directories.examples(),
            Directories.app_annotated_inputs(),
        ):
            for x in directory.glob("*.anno*"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def get_immunopeptidome_options():
        options_dict = {}
        for directory in (
            Directories.immunopeptidomes_humans(),
            Directories.app_immunopeptidomes(),
        ):
            for x in directory.glob("*.bed"):
                options_dict[x.name] = x
        return options_dict

    @staticmethod
    def get_coordinate_options():
        options_dict = {None: None}
        for x in Directories.app_coordinate_files().glob("*.bed"):
            options_dict[x.name] = x
        return options_dict

    @staticmethod
    def get_hla_options():
        hla_types_path = Directories.examples("TCGA_hlaTypesAll.tsv")
        options = pipe(
            ["cut", "-f3", hla_types_path.as_posix()],
            ["tr", ",", "\n"],
            ["sort", "-u"],
        ).split("\n")
        return options

    @staticmethod
    def get_transcript_id_options():
        hla_binders_path = Directories.data(
            "allhlaBinders_exprmean1.IEDBpeps.bed.unique_ids"
        )

        with open(hla_binders_path, "r") as f:
            transcript_options = f.read()

        return transcript_options.split("\n")
