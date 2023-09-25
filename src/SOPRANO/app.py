import os
import pathlib

import streamlit as st

# Find genome options
_HOMO_SAPIENS_DIR = pathlib.Path(__file__).parent / "data" / "homo_sapiens"

_GENOME_DIRS = [item for item in _HOMO_SAPIENS_DIR.glob("*") if item.is_dir()]

_GENOME_NAMES = [
    "{} - Ensembl release {}".format(*x.name.split("_")[::-1])
    for x in _GENOME_DIRS
]

_GENOME_DICT = {
    name: dir_path for name, dir_path in zip(_GENOME_NAMES, _GENOME_DIRS)
}

_ANNO_DIR = pathlib.Path(__file__).parent / "examples"
_ANNO_OPTIONS = {x.name: x for x in _ANNO_DIR.glob("*.anno*")}

_BED_DIR = pathlib.Path(__file__).parent / "immunopeptidomes" / "human"
_BED_OPTIONS = {x.name: x for x in _BED_DIR.glob("*.bed")}

# Init app
st.title("SOPRANO")

# Derived genome definitions
genome_selection = st.selectbox(
    "Select a reference genome:", _GENOME_DICT.keys()
)
ref_id, rel_id = genome_selection.split(" - Ensembl release ")
toplevel_path = (
    _GENOME_DICT[genome_selection] / f"Homo_sapiens.{ref_id}.dna.toplevel.fa"
)
st.text(f"Selected: {toplevel_path}")

# VEP annotated file
input_selection = st.selectbox(
    "Select a VEP annotated file:", _ANNO_OPTIONS.keys()
)
input_path = _ANNO_OPTIONS[input_selection]
st.text(f"Selected: {input_path}")

# BED protein transcript file
bed_selection = st.selectbox("Select a BED protein file:", _BED_OPTIONS.keys())
bed_path = _BED_OPTIONS[bed_selection]
st.text(f"Selected: {bed_path}")

# Substitution method
subs_selection = st.selectbox("Select a substitution method:", (192, 7))
use_ssb192 = subs_selection == 192

# Exclude driver genes
exclude_drives = st.checkbox("Exclude driver genes:", value=True)

# Use randomization
use_random = st.checkbox("Use randomization:", value=False)

if use_random:
    # Select random seed
    random_seed = st.number_input(
        "Select random seed for randomization:", min_value=-1, value="min"
    )

    # TODO: Option for randomization regions

# Pipeline job name & cache
job_name = st.text_input("Define a name for the output of your analysis:")
_DEFAULT_CACHE = pathlib.Path(__file__).parent.parent.parent / "pipeline_cache"
_CACHE = os.environ.get("SOPRANO_CACHE", _DEFAULT_CACHE)
cache_dir = _CACHE / job_name

# If all ok ...
cache_dir.mkdir(exist_ok=True)
