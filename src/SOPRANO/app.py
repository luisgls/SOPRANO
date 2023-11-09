import os
import pathlib

import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from SOPRANO.core import objects
from SOPRANO.utils.app_utils import (
    AnnotatorUIOptions,
    AnnotatorUIProcessing,
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
        genome_ready = genome_selection is not None

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

        cache_selected = st.text_input(
            "Cache directory:", value=Directories.cache().as_posix()
        )
        cache_processed = PipelineUIProcessing.cache(cache_selected)
        cache_ready = os.path.exists(cache_processed)

        job_name_selection = st.text_input(
            "Define a name for the output of your analysis:"
        )
        name_ready = job_name_selection != ""
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

        if st.button(
            "Run Pipeline",
            disabled=not (cache_ready and name_ready and genome_ready),
        ):
            RunTab.pipeline(params=params)


def with_tab_genomes(tab: DeltaGenerator):
    with tab:
        st.subheader("Link existing reference genome files")
        st.markdown(
            "SOPRANO uses a self-contained directory structure for the "
            "the download and caching of derived genomic data.\n\n"
            "Users who have an existing Ensembl VEP configuration on their "
            "computer may have readily available reference genomes "
            "downloaded. These downloads can be linked to the SOPRANO "
            "directories by running the button below.\n\nNote: "
            " the default VEP cache that is searched for is `~/.vep` but you "
            "can define any other non-standard locations that reference "
            "genomes may be found within."
        )

        cache_location_selection = st.text_input(
            "VEP cache location:", value=Directories.std_sys_vep().as_posix()
        )

        cache_location_processed = LinkVEPUIProcessing.cache_location(
            cache_location_selection
        )

        if st.button("Attempt VEP cache link", disabled=False):
            RunTab.link_vep(cache_location_processed)

        st.subheader("Download new reference genome files")
        st.markdown(
            "You can use SOPRANO to download reference genomes from the "
            "ensembl FTP server by making use of the below definitions, before"
            " clicking `Download`.\n\n"
            "The core SOPRANO calculation requires toplevel reference data "
            "to be available. Secondary downloads of the primary assembly "
            "may be used to accelerate the annotation procedure; though this "
            "is _not_ essential."
        )

        species_selection = st.text_input(
            "Define the species",
            value="Homo Sapiens",
            disabled=True,
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
        st.markdown(
            "Generate an annotated mutation file from a directory containing "
            "VCF files suitable for consumption by SOPRANO."
        )

        vcf_dir_selection = st.text_input(
            "Directory containing VCF files:", value=pathlib.Path.home()
        )
        vcf_dir_ready, vcf_dir_processed = AnnotatorUIProcessing.vcf_dir(
            vcf_dir_selection
        )

        assembly_selection = st.selectbox(
            "Genome reference assembly:",
            options=AnnotatorUIOptions.genome_assembly(),
        )
        (
            assembly_ready,
            assembly_processed,
        ) = AnnotatorUIProcessing.genome_assembly(assembly_selection)

        name_selection = st.text_input(
            "Choose a name for the annotated output:"
        )
        name_ready, name = AnnotatorUIProcessing.output_name(name_selection)

        if st.button(
            "Annotate",
            disabled=not (vcf_dir_ready and assembly_ready and name_ready),
        ):
            RunTab.annotate(sources_dir=vcf_dir_processed, output_name=name)


def with_tab_info(tab: DeltaGenerator):
    with tab:
        st.title("Welcome to SOPRANO! :wave:")
        st.caption("Selection On PRotein ANnotated regiOns")
        st.markdown(
            "This application is designed to provide a user interface to the "
            "SOPRANO computational pipeline, without the need of command line "
            "intervention."
            "\n\n"
            "There are three essential files required to run "
            "SOPRANO. These define the\n"
            "1. Reference genome\n"
            "2. Annotated somatic mutations\n"
            "3. Immunopeptidome\n"
            "\n\n"
            "These three inputs can be configured in term via the tabs "
            "indicating steps 1, 2 and 3. Once you have prepared your data, "
            "step 4 will enable you to run the pipeline, subject to "
            "further runtime configuration choices."
            "\n\n"
            "Any technical issues can be raised on [GitHub]"
            "(https://github.com/instituteofcancerresearch/SOPRANO/issues)"
        )


def with_tab_immunopeptidome(tab: DeltaGenerator):
    with tab:
        st.header("HLA-allele selections")

        st.markdown(
            "In order to generate custom, or, patient specific "
            "immunopeptidome files to run through SOPRANO, we must "
            "specify a set of HLA alleles to restrict the analysis to.\n\n"
            "You should enter 1-6 HLA allele choices into the text box,"
            " or upload a text file containing similar information."
        )

        hla_ready, alleles_selected = text_or_file(
            "Define the HLA allele selection manually via text input *OR* "
            "file uploader. (1-6 alleles are required).",
            min_args=1,
            max_args=6,
            help_text="Provide HLA alleles on separate lines.",
            help_upload="Select a text file defining the HLA alleles on "
            "separate lines.",
        )

        if not hla_ready:
            st.warning("HLA selection is currently invalid.")

        alleles_processed = ImmunopeptidomeUIProcessing.hla_alleles(
            [] if alleles_selected is None else alleles_selected
        )

        st.header("Ensembl transcript selections")

        st.markdown(
            "Analyses can be further restricted by providing a set of "
            "Ensembl transcript IDs to either\n\n"
            "1. Exclude from the analysis; or\n"
            "2. Restrict the analysis to.\n\n"
            "These filtering characteristics are presented as the choices\n"
            "- 'Retention'\n"
            "- 'Exclusion'\n\n"
            "in the following drop down menu.\n\n"
            "Transcript IDs should either be manually entered into the text "
            "box over separate lines, or uploaded within a text file of "
            "similar contents."
        )

        transcripts_ready, transcripts_selected = text_or_file(
            "Define the Ensembl transcript IDs to restrict or exclude from the"
            " immunopeptidome construction.",
            help_text="Provide transcript IDs on separate lines.",
            help_upload="Select a text file defining transcript IDs on "
            "separate lines.",
        )

        transcripts_processed = ImmunopeptidomeUIProcessing.transcript_ids(
            [] if transcripts_selected is None else transcripts_selected
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

        st.markdown(
            "Once you are happy with your immunopeptidome choices, "
            "provide a name for the corresponding file, then click "
            "'Build immunopeptidome'."
        )

        name_selected = st.text_input("Immunopeptidome name:")
        name_processed = ImmunopeptidomeUIProcessing.name(name_selected)
        name_ready = name_processed != ".bed"

        if not name_ready:
            st.warning("Please provide an input file name.")

        if st.button(
            "Build immunopeptidome",
            disabled=not (hla_ready and name_ready),
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
            "Step 1: Prepare genome reference",
            "Step 2: Annotate mutations",
            "Step 3: Prepare immunopeptidome",
            "Step 4: Run pipeline",
        ]
    )
    with_tab_info(welcome_tab)
    with_tab_genomes(genome_tab)
    with_tab_annotator(annotate_tab)
    with_tab_immunopeptidome(immunopeptidome_tab)
    with_tab_pipeline(pipeline_tab)
