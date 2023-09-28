import argparse
import pathlib
import time

from shiny import App, reactive, render, ui

from SOPRANO import objects, run_local_ssb_selection
from SOPRANO.misc_utils import Directories

_CACHE = Directories.cache()

# Find genome options
_HOMO_SAPIENS_DIR = Directories.homo_sapien_genomes()

_GENOME_DIRS = [item for item in _HOMO_SAPIENS_DIR.glob("*") if item.is_dir()]

_GENOME_NAMES = [
    "{} - Ensembl release {}".format(*x.name.split("_")[::-1])
    for x in _GENOME_DIRS
]

_GENOME_DICT = {
    name: dir_path for name, dir_path in zip(_GENOME_NAMES, _GENOME_DIRS)
}

_ANNO_DIR = Directories.examples()
_ANNO_OPTIONS = {x.name: x for x in _ANNO_DIR.glob("*.anno*")}

_BED_DIR = Directories.immuno_humans()
_BED_OPTIONS = {x.name: x for x in _BED_DIR.glob("*.bed")}

_SUB_OPTIONS = ("SSB7", "SSB192")


def process_genome_selection(genome_selection: str) -> tuple:
    ref_id, rel_id = genome_selection.split(" - Ensembl release ")
    toplevel_path = (
        _GENOME_DICT[genome_selection]
        / f"Homo_sapiens.{ref_id}.dna.toplevel.fa"
    )
    return toplevel_path, ref_id, rel_id


def process_input_selection(input_selection: str) -> pathlib.Path:
    return _ANNO_OPTIONS[input_selection]


def process_bed_selection(bed_selection: str) -> pathlib.Path:
    return _BED_OPTIONS[bed_selection]


def process_subs_selection(subs_selection: str) -> bool:
    return subs_selection == "SSB192"


def process_job_name(job_name: str) -> pathlib.Path:
    return _CACHE.joinpath(job_name)


def _options_ui(opt_id, opt_desc, opt_vals, out_id):
    return (
        ui.input_select(opt_id, opt_desc, opt_vals),
        ui.output_text_verbatim(out_id),
    )


app_ui = ui.page_fluid(
    ui.panel_title("SOPRANO"),
    *_options_ui(
        "genome_selection",
        "Reference Genome",
        tuple(_GENOME_DICT.keys()),
        "g_selection",
    ),
    *_options_ui(
        "input_selection",
        "VEP Annotated File",
        tuple(_ANNO_OPTIONS.keys()),
        "i_selection",
    ),
    *_options_ui(
        "bed_selection",
        "BED Protein File",
        tuple(_BED_OPTIONS.keys()),
        "b_selection",
    ),
    *_options_ui(
        "subs_selection", "Substitution Method", _SUB_OPTIONS, "s_selection"
    ),
    ui.input_checkbox("exclude_drivers", "Exclude Driver Genes", value=True),
    ui.input_checkbox("use_random", "Use Randomization", value=False),
    ui.input_numeric("random_seed", "Random Seed", value=-1, min=-1),
    ui.input_text("job_name", "Analysis Name"),
    ui.output_text_verbatim("c_location"),
    ui.input_action_button("run_pipeline", "Run Pipeline"),
    ui.output_text_verbatim("result", placeholder=True),
)


def server(input, output, session):
    @output
    @render.text
    def g_selection():
        return process_genome_selection(input.genome_selection())[1]

    @output
    @render.text
    def i_selection():
        return process_input_selection(input.input_selection())

    @output
    @render.text
    def b_selection():
        return process_bed_selection(input.bed_selection())

    @output
    @render.text
    def c_location():
        return process_job_name(input.job_name())

    @output
    @render.text
    @reactive.event(input.run_pipeline)
    def result():
        namespace = argparse.Namespace(
            analysis_name=input.job_name(),
            input_path=process_input_selection(input.input_selection()),
            bed_path=process_bed_selection(input.bed_selection()),
            cache_dir=process_job_name(input.job_name()),
            target_regions=None,  # TODO: Fix
            use_ssb192=process_subs_selection(input.subs_selection()),
            use_random=input.use_random(),
            exclude_drivers=input.exclude_drivers(),
            seed=input.random_seed(),
            transcript=objects.EnsemblTranscripts.transcript_length,
            protein_transcript=objects.EnsemblTranscripts.protein_transcript_length,
            transcript_ids=objects.EnsemblTranscripts.transcript_fasta,
            genome_ref=process_genome_selection(input.genome_selection())[1],
            release=process_genome_selection(input.genome_selection())[2],
        )

        params = objects.Parameters.from_namespace(namespace)
        params.cache_dir.mkdir(exist_ok=True)

        t_start = time.time()
        run_local_ssb_selection.main(namespace)
        t_end = time.time()

        return f"Pipeline run in {int(t_end - t_start)} seconds"


app = App(app_ui, server)
