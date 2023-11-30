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
  file.path(translator_dir, "ENSP2ENST.txt"), delim = "\t", col_names = TRUE, show_col_types = FALSE
)
ref2vep_data <- readr::read_delim(
  file.path(translator_dir, "REF2VEP.txt"), delim = "\t", col_names = TRUE, show_col_types = FALSE
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

  print(bialleilic_indels_extracted)

  gt_tidied <- extract_gt_tidy(bialleilic_indels_extracted)

  print(gt_tidied)

  fields_fixed <- fix_vcf_fields(bialleilic_indels_extracted)

  print(fields_fixed)

  stop()

  merged <- merge_gt_tidy_with_fixed_fields(gt_tidied, fields_fixed)
  reduced <- reduce_tibble_id_to_chr_pos_ref_alt(merged)
}

for (path in vcf_file_paths) {
  annotate_vcf_file(
    path, NULL, ensp2enst_data, ref2vep_data,
    cv = covs, refdb = refdb
  )
}


# print(vcf_file_paths)
# stop()
#
#
# # Read mutations from VCF files
# vcf_names <- c()
# for (path in vcf_files) {
#   vcf_names <- c(vcf_names, basename(path))
# }
#
#
# # vcf_data <- lapply(vcf_files, function(vcf) {
# #   v <- vcfR::read.vcfR(vcf)
# #   v <- vcfR::extract.indels(v)
# #   v <- v[vcfR::is.biallelic(v), ]
# # })
# # names(vcf_data) <- vcf_names
# #
# # gt_fields <- lapply(vcf_data, function(vcf) {
# #   vcfR::extract_gt_tidy(vcf)
# # })
# #
# # fix_fields <- lapply(vcf_data, function(vcf) {
# #   vcf@fix %>% tibble::as_tibble()
# # })
# #
# # merged_list <- mapply(c, fix_fields, gt_fields, SIMPLIFY = FALSE)
# #
# # # Create table with basic data from VCF
# # all_with_vaf <- lapply(merged_list, function(test) {
# #   test2 <- tibble::as_tibble(cbind(test$CHROM, test$POS, test$REF, test$ALT))
# # })
#
# # df_all_with_vaf <- dplyr::bind_rows(all_with_vaf, .id = "name")
# # names(df_all_with_vaf) <- c(
# #   "sampleID", "chr", "position", "ref_allele", "alt_allele"
# # )
#
# # Read transcript and annotation info
# transcriptlist <- readr::read_delim(
#   ensp_2_enst_path, delim = "\t", col_names = TRUE
# )
# variantlist <- readr::read_delim(
#   variant_list_path, delim = "\t", col_names = TRUE
# )
#
# # See discussion:
# # https://github.com/im3sanger/dndscv/issues/30#issuecomment-1000868593
# #
# # load(covs_path) # Loads the covs object
# #
# #
# # # dndscv bit...
# #
# # # Run dndscv with GRCh38 defs... required a bit of tinkering
# # # Other genome need to download Rda file-> see dndscv tutorial website.
# #
# # # See discussion:
# # # unix.stackexchange.com/questions/497990/
# # #   how-to-know-if-rsync-did-not-change-any-files
# # df_all_with_vaf$chr <- gsub("chr", "", as.vector(df_all_with_vaf$chr))
# #
# # if (assembly == "GRCh38") {
# #   res1dnds <- dndscv::dndscv(
# #     mutations = df_all_with_vaf,
# #     outmats = TRUE,
# #     max_muts_per_gene_per_sample = Inf,
# #     max_coding_muts_per_sample = Inf,
# #     outp = 2,
# #     use_indel_sites = TRUE,
# #     min_indels = 1,
# #     refdb = refdb_path,
# #     cv = covs
# #   )
# # } else {
# #   res1dnds <- dndscv::dndscv(
# #     mutations = df_all_with_vaf,
# #     outmats = TRUE,
# #     max_muts_per_gene_per_sample = Inf,
# #     max_coding_muts_per_sample = Inf,
# #     outp = 2,
# #     use_indel_sites = TRUE,
# #     min_indels = 1
# #   )
# # }
# #
# #
# # ## Get table for annotated mutations
# # annotation <- dplyr::as_tibble(res1dnds$annotmuts)
# #
# # ## Join protein ID with transcript ID
# # annot1 <- annotation %>%
# #   dplyr::left_join(., transcriptlist, by = c("pid" = "ProteinstableID")) %>%
# #   dplyr::left_join(., variantlist, by = c("impact" = "REF")) %>%
# #   dplyr::select(
# #     chr, pos, ref, mut, TranscriptstableID, VEP, aachange, ntchange, codonsub
# #   )
# #
# # ## Parse data
# # annot1$codonsub2 <- gsub(">", "/", annot1$codonsub)
# # annot1$aachange2 <- gsub("[0-9]+", "/", annot1$aachange)
# #
# # number <- "[0-9]+"
# #
# # annot2 <- annot1 %>%
# #   tidyr::unite(idtmp, c("chr", "pos"), sep = "_") %>%
# #   tidyr::unite(change, c("ref", "mut"), sep = "/") %>%
# #   tidyr::unite(id, c("idtmp", "change"), sep = "_") %>%
# #   dplyr::mutate(
# #     protpos = stringr::str_extract(annot1$aachange, number),
# #     cdspos = stringr::str_extract(annot1$ntchange, number)
# #   )
# #
# # dummy_ssb <- annot2 %>%
# #   dplyr::mutate(
# #     col1 = "NA", col2 = "NA", col3 = "NA", col4 = "Transcript", col5 = "NA",
# #     col6 = "NA", col7 = "NA", col8 = "NA"
# #   ) %>%
# #   dplyr::select(
# #     id, col1, col2, col3, TranscriptstableID,
# #     col4, VEP, col5, cdspos, protpos, aachange2, codonsub2
# #   )
# #
# # dummy_ssb %<>% dplyr::mutate(aachange2 = as.character(aachange2))
# # dummy_ssb %<>% dplyr::mutate(codonsub2 = as.character(codonsub2))
# #
# # df_ssb <- as.data.frame(dummy_ssb)
# #
# # df_ssb <- df_ssb[!is.na(df_ssb$TranscriptstableID),]
# #
# # write.table <(
# #   df_ssb,
# #   file = out_name, quote = FALSE, sep = "\t",
# #   col.names = FALSE, row.names = FALSE
# # )
