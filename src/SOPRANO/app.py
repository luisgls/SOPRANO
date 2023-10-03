import time

import pandas as pd
import streamlit as st

import SOPRANO.utils.path_utils
from SOPRANO.core import objects
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.app_utils import get_human_genome_options, st_capture
from SOPRANO.utils.path_utils import Directories

genome_options = get_human_genome_options()

_ANNO_DIR = Directories.examples()
_ANNO_OPTIONS = {x.name: x for x in _ANNO_DIR.glob("*.anno*")}

_BED_DIR = Directories.immuno_humans()
_BED_OPTIONS = {x.name: x for x in _BED_DIR.glob("*.bed")}


if __name__ == "__main__":
    # Init app
    st.title("SOPRANO")
    st.caption("Selection On PRotein ANnotated regiOns")

    def process_genome_selection():
        genome_selection = st.session_state.genome_selection

        ref_id, rel_id = genome_selection.split(" - Ensembl release ")

        (
            genomes_path,
            chroms_path,
        ) = SOPRANO.utils.path_utils.genome_pars_to_paths(ref_id, rel_id)

        st.session_state.ref_genome = objects.GenomePaths(
            sizes=chroms_path, fasta=genomes_path
        )

        st.text(f"Selected: {st.session_state.ref_genome.fasta}")

    # Derived genome definitions
    st.selectbox(
        "Select a reference genome:",
        genome_options.keys(),
        key="genome_selection",
    )
    process_genome_selection()

    def process_annotation():
        input_selection = st.session_state.input_selection
        st.session_state.input_path = _ANNO_OPTIONS[input_selection]
        st.text(f"Selected: {st.session_state.input_path}")

    # VEP annotated file
    st.selectbox(
        "Select a VEP annotated file:",
        _ANNO_OPTIONS.keys(),
        key="input_selection",
    )
    process_annotation()

    def process_bed():
        bed_selection = st.session_state.bed_selection
        st.session_state.bed_path = _BED_OPTIONS[bed_selection]
        st.text(f"Selected: {st.session_state.bed_path}")

    # BED protein transcript file
    st.selectbox(
        "Select a BED protein file:", _BED_OPTIONS.keys(), key="bed_selection"
    )
    process_bed()

    def process_subs():
        subs_selection = st.session_state.subs_selection
        st.session_state.use_ssb192 = subs_selection == 192
        st.text(f"Using SSB192: {st.session_state.use_ssb192}")

    # Substitution method
    st.selectbox(
        "Select a substitution method:", (192, 7), key="subs_selection"
    )
    process_subs()

    # Exclude driver genes
    st.checkbox("Exclude driver genes:", value=True, key="exclude_drivers")

    # Use randomization
    st.checkbox("Use randomization:", value=False, key="use_random")

    if st.session_state.use_random:
        # Select random seed
        st.number_input(
            "Select random seed for randomization:",
            min_value=-1,
            value="min",
            key="random_seed",
        )
        # TODO: Option for randomization regions
    else:
        st.session_state.random_seed = -1

    def process_name():
        job_name = st.session_state.job_name
        st.session_state.cache_dir = Directories.cache(job_name)
        st.text(f"Cache location: {st.session_state.cache_dir}")

    # Pipeline job name & cache
    st.text_input(
        "Define a name for the output of your analysis:", key="job_name"
    )
    process_name()

    st.session_state.params = objects.Parameters(
        analysis_name=st.session_state.job_name,
        input_path=st.session_state.input_path,
        bed_path=st.session_state.bed_path,
        cache_dir=st.session_state.cache_dir,
        random_regions=None,  # TODO: Fix
        use_ssb192=st.session_state.use_ssb192,
        use_random=st.session_state.use_random,
        exclude_drivers=st.session_state.exclude_drivers,
        seed=st.session_state.random_seed,
        transcripts=objects.TranscriptPaths.defaults(),
        genomes=st.session_state.ref_genome,
    )

    st.session_state.job_complete = False

    def run_pipeline_in_app():
        st.session_state.cache_dir.mkdir(exist_ok=True)

        t_start = time.time()
        output = st.empty()
        with st_capture(output.code):
            run_pipeline(st.session_state.params)
            # run_local_ssb_selection.main(st.session_state.namespace)
        t_end = time.time()

        st.session_state.compute_time_str = (
            f"Pipeline run in {int(t_end - t_start)} seconds"
        )

        st.session_state.data_frame = pd.read_csv(
            st.session_state.params.results_path, sep="\t"
        )

        st.session_state.job_complete = True

        st.text(st.session_state.compute_time_str)
        st.dataframe(st.session_state.data_frame, hide_index=True)
        st.text(f"dN/dS file: {st.session_state.params.results_path}")

    if st.button("Run Pipeline"):
        run_pipeline_in_app()
