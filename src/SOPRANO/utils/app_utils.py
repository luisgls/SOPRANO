from contextlib import contextmanager, redirect_stdout
from io import StringIO

from SOPRANO.utils.path_utils import Directories


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


def get_human_genome_options():
    homo_sapiens_dir = Directories.homo_sapien_genomes()

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
