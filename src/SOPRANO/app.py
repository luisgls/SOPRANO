import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from SOPRANO.core import objects
from SOPRANO.utils.app_utils import (
    DownloaderUIOptions,
    DownloaderUIProcessing,
    ImmunopeptidomesUIOptions,
    ImmunopeptidomeUIProcessing,
    LinkVEPUIProcessing,
    PipelineUIOptions,
    PipelineUIProcessing,
    RunTab,
    text_or_file,
)
from SOPRANO.utils.path_utils import Directories


def with_tab_pipeline(tab: DeltaGenerator):
    with tab:
        st.title("SOPRANO")
        st.caption("Selection On PRotein ANnotated regiOns")

        genome_selection = st.selectbox(
            "Select a reference genome:",
            PipelineUIOptions.genome_reference(),
        )
        genome_processed = PipelineUIProcessing.genome_reference(
            genome_selection
        )

        annotation_selection = st.selectbox(
            "Select a VEP annotated file:",
            PipelineUIOptions.annotated_mutations(),
        )
        annotation_processed = PipelineUIProcessing.annotated_mutations(
            annotation_selection
        )

        immunopeptidome_selection = st.selectbox(
            "Select a BED protein file:",
            PipelineUIOptions.immunopeptidome(),
        )
        immunopeptidome_processed = PipelineUIProcessing.immunopeptidome(
            immunopeptidome_selection
        )

        substitution_selection = st.selectbox(
            "Select a substitution method:",
            PipelineUIOptions.substitution_method(),
        )
        substitution_processed = PipelineUIProcessing.substitution_method(
            substitution_selection
        )

        exclude_drivers = st.checkbox("Exclude driver genes:", value=True)
        use_randomization = st.checkbox("Use randomization:", value=False)

        if use_randomization:
            random_seed = st.number_input(
                "Select random seed for randomization:",
                min_value=-1,
                value="min",
            )
            coordinates_selection = st.selectbox(
                "Select a BED file defining the regions to randomize over:",
                PipelineUIOptions.coordinates(),
            )
            coordinates_processed = PipelineUIProcessing.coordinates(
                coordinates_selection
            )
        else:
            random_seed = -1
            coordinates_processed = None

        job_name_selection = st.text_input(
            "Define a name for the output of your analysis:"
        )
        job_name_processed = PipelineUIProcessing.job_name(job_name_selection)

        params = objects.Parameters(
            analysis_name=job_name_selection,
            input_path=annotation_processed,
            bed_path=immunopeptidome_processed,
            cache_dir=job_name_processed,
            random_regions=coordinates_processed,
            use_ssb192=substitution_processed == 192,
            use_random=use_randomization,
            exclude_drivers=exclude_drivers,
            seed=random_seed,
            transcripts=objects.TranscriptPaths.defaults(),
            genomes=genome_processed,
        )

        if st.button("Run Pipeline"):
            RunTab.pipeline(params=params)


def with_tab_genomes(tab: DeltaGenerator):
    with tab:
        st.title("Download Reference Genomes")

        st.caption(
            f"Symbolically link files in your VEP cache to the SOPRANO data "
            f"folder: {Directories.data()}"
        )

        cache_location_selection = st.text_input(
            "VEP cache location:", value=Directories.std_sys_vep().as_posix()
        )
        cache_location_processed = LinkVEPUIProcessing.cache_location(
            cache_location_selection
        )

        if st.button("Link", disabled=False):
            RunTab.link_vep(cache_location_processed)

        st.text(
            "Download Ensembl genome reference data. See "
            "https://www.ensembl.org/info/genome/variation/species/"
            "species_data_types.html for definitions."
        )

        species_selection = st.text_input(
            "Define the species",
            value="Homo Sapiens",
        )
        species_processed = DownloaderUIProcessing.species(species_selection)

        assembly_selection = st.text_input(
            "Define the genome reference:",
            value="GRCh38",
        )
        assembly_processed = DownloaderUIProcessing.assembly(
            assembly_selection
        )

        release_selection = st.number_input(
            "Define the Ensembl release:",
            min_value=76,
            key="download_release",
            value=110,
        )
        release_processed = DownloaderUIProcessing.release(release_selection)

        type_selection = st.selectbox(
            "Download type:",
            options=DownloaderUIOptions.type(),
        )
        type_processed = DownloaderUIProcessing.type(type_selection)

        if st.button("Download", disabled=False):
            RunTab.download(
                species=species_processed,
                assembly=assembly_processed,
                release=release_processed,
                download_type=type_processed,
            )


def with_tab_annotator(tab: DeltaGenerator):
    with tab:
        st.title("Annotate VCF File")
        st.caption(
            "Upload a VCF file to annotate for use in the SOPRANO pipeline."
        )
        st.file_uploader("Select a VCF file:", key="vcf_selection")

        if st.button("Annotate", disabled=True):
            RunTab.annotate()


def with_tab_info(tab: DeltaGenerator):
    with tab:
        st.title("Information")
        st.caption("Description of what is going on...")


def with_tab_immunopeptidome(tab: DeltaGenerator):
    with tab:
        st.title("Immunopeptidomes")
        st.caption(
            "Generate a patient specific immunopeptidome file derived from "
            "the SOPRANO internal master file and HLA allele selections."
        )

        ready, other = text_or_file("This is for IP")

        st.text(f"Status: {ready}")
        st.text(f"{other}, {type(other)}")

        alleles_selected = st.multiselect(
            "Select HLA alleles (min 2, max 6):",
            options=ImmunopeptidomesUIOptions.hla_alleles(),
            max_selections=6,
        )
        alleles_processed = ImmunopeptidomeUIProcessing.hla_alleles(
            alleles_selected
        )

        name_selected = st.text_input("Immunopeptidome name:")
        name_processed = ImmunopeptidomeUIProcessing.name(name_selected)

        transcripts_selected = st.multiselect(
            "Select transcript IDs for subset (optional):",
            options=ImmunopeptidomesUIOptions.transcript_ids(),
        )
        transcripts_processed = ImmunopeptidomeUIProcessing.transcript_ids(
            transcripts_selected
        )

        subset_method_selected = st.selectbox(
            "Select method to subset available Ensembl transcript IDs "
            "(optional):",
            options=ImmunopeptidomesUIOptions.subset_method(),
        )
        (
            transcripts_retained,
            transcripts_excluded,
        ) = ImmunopeptidomeUIProcessing.subset_method(
            transcripts_processed, subset_method_selected
        )

        if st.button(
            "Build immunopeptidome",
            disabled=not (2 <= len(alleles_processed) <= 6),
        ):
            RunTab.immunopeptidome(
                alleles_processed,
                output_name=name_processed,
                transcripts_retained=transcripts_retained,
                transcripts_excluded=transcripts_excluded,
            )


if __name__ == "__main__":
    st.set_page_config(layout="wide")
    (
        welcome_tab,
        genome_tab,
        annotate_tab,
        immunopeptidome_tab,
        pipeline_tab,
    ) = st.tabs(
        [
            "Welcome!",
            "Step 1: Prepare genome references",
            "Step 2: Annotate mutation files",
            "Step 3: Prepare immunopeptidome files",
            "Step 4: Run pipeline",
        ]
    )
    with_tab_info(welcome_tab)
    with_tab_genomes(genome_tab)
    with_tab_annotator(annotate_tab)
    with_tab_immunopeptidome(immunopeptidome_tab)
    with_tab_pipeline(pipeline_tab)
