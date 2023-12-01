library(magrittr)

source("check_args.R")
source("vcf_gymnastics.R")

# Define parser
opt_parser <- optparse::OptionParser()

# Data directory containing vcf.gz files
opt_parser <- optparse::add_option(
  opt_parser, c("-s", "--source"),
  type = "character",
  help = "Path to .vcf.gz file, or directory containing .vcf.gz files.",
  metavar = "character",
)

# Output filename
opt_parser <- optparse::add_option(
  opt_parser, c("-o", "--output"),
  type = "character", default = NULL,
  help = "Output file path.", metavar = "character"
)

# Location of sources for translating ensp and ref
opt_parser <- optparse::add_option(
  opt_parser, c("-t", "--translate"),
  type = "character",
  help = "Location of translator files.", metavar = "character"
)

# Genome assembly
opt_parser <- optparse::add_option(
  opt_parser, c("-a", "--assembly"),
  type = "character", default = "GRCh38",
  help = "Homo sapiens genome assembly", metavar = "character"
)

# Parse inputs
args <- optparse::parse_args(opt_parser)
source_path <- args$source
output_path <- args$output
translator_dir <- args$translate
assembly_type <- args$assembly

# Check inputs
check_source_path(source_path)
check_output_path(output_path)
check_translator_dir(translator_dir)
check_assembly_type(assembly_type)
check_auxiliary_paths(
  translator_dir,
  "ENSP2ENST.txt",
  "REF2VEP.txt",
  "covariates_hg19_hg38_epigenome_pcawg.rda",
  "RefCDS_human_GRCh38_GencodeV18_recommended.rda"
)

# Load auxiliary data
ensp2enst_data <- readr::read_delim(
  file.path(translator_dir, "ENSP2ENST.txt"),
  delim = "\t", col_names = TRUE, show_col_types = FALSE
)
ref2vep_data <- readr::read_delim(
  file.path(translator_dir, "REF2VEP.txt"),
  delim = "\t", col_names = TRUE, show_col_types = FALSE
)
if (assembly_type == "GRCh38") {
  # sets covs object
  load(file.path(translator_dir, "covariates_hg19_hg38_epigenome_pcawg.rda"))
  refdb <- file.path(
    translator_dir, "RefCDS_human_GRCh38_GencodeV18_recommended.rda"
  )
} else {
  covs <- "hg19"
  refdb <- "hg19"
}

# Select vcf.gz files from input source
vcf_file_paths <- select_vcf_files(source_path)

annotate_vcf_file <- function(vcf_file_path, output_path, ensp2enst, ref2vep, cv, refdb) {
  if (is.null(output_path)) {
    dir_name <- dirname(vcf_file_path)
    file_name <- paste0(tools::file_path_sans_ext(basename(vcf_file_path)), ".anno")
    output_path <- file.path(dir_name, file_name)
  }

  bialleilic_indels_extracted <- extract_biallelic_indels(vcf_file_path)

  gt_tidied <- extract_gt_tidy(bialleilic_indels_extracted)

  fields_fixed <- fix_vcf_fields(bialleilic_indels_extracted)

  merged <- merge_gt_tidy_with_fixed_fields(gt_tidied, fields_fixed)

  reduced <- reduce_tibble_id_to_chr_pos_ref_alt(merged)

  chr_patterns_corrected <- correct_chr_patterns(reduced)

  sample_name <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(vcf_file_path)))

  chr_patterns_corrected$sampleID <- sample_name

  dndscv_reuslt <- dndscv::dndscv(
    mutations = chr_patterns_corrected,
    outmats = TRUE,
    max_muts_per_gene_per_sample = Inf,
    max_coding_muts_per_sample = Inf,
    outp = 2,
    use_indel_sites = TRUE,
    min_indels = 1,
    refdb = refdb,
    cv = covs
  )

  with_transcript_and_variant_info <- add_transcript_and_variant_info(dndscv_reuslt, transcriptlist = ensp2enst, variantlist = ref2vep)

  with_reformatted_annotations <- update_annotation_formatting(with_transcript_and_variant_info)

  prepare_ssb_input_file(with_reformatted_annotations, output_path)
}

for (path in vcf_file_paths) {
  annotate_vcf_file(
    path, NULL, ensp2enst_data, ref2vep_data,
    cv = covs, refdb = refdb
  )
}
