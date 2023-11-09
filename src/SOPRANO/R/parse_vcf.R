library(magrittr)

# Define parser
opt_parser <- optparse::OptionParser()

# Data directory containing vcf.gz files
opt_parser <- optparse::add_option(
  opt_parser, c("-d", "--dir"),
  type = "character",
  help = "Directory containing .vcf.gz files.", metavar = "character",
)

# Output filename
opt_parser <- optparse::add_option(
  opt_parser, c("-o", "--out"),
  type = "character", default = NULL,
  help = "Output file path.", metavar = "character"
)

# Location of sources for translating ensp and ref
opt_parser <- optparse::add_option(
  opt_parser, c("-t", "--trans"),
  type = "character",
  help = "Location of translator files.", metavar = "character"
)

# Parse inputs
args <- optparse::parse_args(opt_parser)

# Get data sources
data_dir <- args$dir

if (is.null(data_dir)) {
    stop("Data directory not defined. Flag -d | --d <directory path> required")
}
if (!dir.exists(data_dir)) {
    stop(paste("Data directory does not exist:", data_dir))
}


# Get translator dir
trans_dir <- args$trans

if (is.null(trans_dir)) {
    stop("Translator directory not defined. Flag -t | --t <dir path> required")
}

if (!dir.exists(trans_dir)) {
    stop(paste("Translator directory does not exits:", trans_dir))
}

ensp_2_enst_path <- file.path(trans_dir, "ENSP2ENST.txt")
variant_list_path <- file.path(trans_dir, "REF2VEP.txt")
covs_path = file.path(trans_dir, "covariates_hg19_hg38_epigenome_pcawg.rda")
refdb_path <- file.path(
    trans_dir, "RefCDS_human_GRCh38_GencodeV18_recommended.rda"
)

for (p in c(ensp_2_enst_path, variant_list_path, covs_path, refdb_path)) {
    if (! file.exists(p)) {
        stop(paste("Auxiliary file not found:", p))
    }
}

# Get output file name
if (is.null(args$out)) {
  out_name <- file.path(data_dir, paste0(basename(data_dir), ".anno"))
} else {
  out_name <- args$out
}

# Read mutations from VCF files
vcf_files <- Sys.glob(file.path(args$dir, "*.vcf.gz"))
vcf_names <- sub(
  "_TAIL_somatic_snvs_snpEff.ann.vcf.gz", "", sub("\\.\\/", "", vcf_files)
)
vcf_data <- lapply(vcf_files, function(vcf) {
  v <- vcfR::read.vcfR(vcf)
  v <- vcfR::extract.indels(v)
  v <- v[vcfR::is.biallelic(v), ]
})
names(vcf_data) <- vcf_names

gt_fields <- lapply(vcf_data, function(vcf) {
  vcfR::extract_gt_tidy(vcf)
})

fix_fields <- lapply(vcf_data, function(vcf) {
  vcf@fix %>% tibble::as_tibble()
})

merged_list <- mapply(c, fix_fields, gt_fields, SIMPLIFY = FALSE)

# Create table with basic data from VCF
all_with_vaf <- lapply(merged_list, function(test) {
  test2 <- tibble::as_tibble(cbind(test$CHROM, test$POS, test$REF, test$ALT))
})

df_all_with_vaf <- dplyr::bind_rows(all_with_vaf, .id = "name")
names(df_all_with_vaf) <- c(
  "sampleID", "chr", "position", "ref_allele", "alt_allele"
)

# Read transcript and annotation info



transcriptlist <- readr::read_delim(
  ensp_2_enst_path, delim = "\t", col_names = TRUE
)
variantlist <- readr::read_delim(
  variant_list_path, delim = "\t", col_names = TRUE
)

# See discussion:
# https://github.com/im3sanger/dndscv/issues/30#issuecomment-1000868593

load(covs_path) # Loads the covs object


# dndscv bit...

# Run dndscv with GRCh38 defs... required a bit of tinkering
# Other genome need to download Rda file-> see dndscv tutorial website.

# See discussion:
# unix.stackexchange.com/questions/497990/
#   how-to-know-if-rsync-did-not-change-any-files
df_all_with_vaf$chr <- gsub("chr", "", as.vector(df_all_with_vaf$chr))


res1dnds <- dndscv::dndscv(
  mutations = df_all_with_vaf,
  outmats = TRUE,
  max_muts_per_gene_per_sample = Inf,
  max_coding_muts_per_sample = Inf,
  outp = 2,
  use_indel_sites = TRUE,
  min_indels = 1,
  refdb = refdb_path,
  cv = covs
)


## Get table for annotated mutations
annotation <- dplyr::as_tibble(res1dnds$annotmuts)

## Join protein ID with transcript ID
annot1 <- annotation %>%
  dplyr::left_join(., transcriptlist, by = c("pid" = "ProteinstableID")) %>%
  dplyr::left_join(., variantlist, by = c("impact" = "REF")) %>%
  dplyr::select(
    chr, pos, ref, mut, TranscriptstableID, VEP, aachange, ntchange, codonsub
  )

## Parse data
annot1$codonsub2 <- gsub(">", "/", annot1$codonsub)
annot1$aachange2 <- gsub("[0-9]+", "/", annot1$aachange)

number <- "[0-9]+"

annot2 <- annot1 %>%
  tidyr::unite(idtmp, c("chr", "pos"), sep = "_") %>%
  tidyr::unite(change, c("ref", "mut"), sep = "/") %>%
  tidyr::unite(id, c("idtmp", "change"), sep = "_") %>%
  dplyr::mutate(
    protpos = stringr::str_extract(annot1$aachange, number),
    cdspos = stringr::str_extract(annot1$ntchange, number)
  )

dummy_ssb <- annot2 %>%
  dplyr::mutate(
    col1 = "NA", col2 = "NA", col3 = "NA", col4 = "Transcript", col5 = "NA",
    col6 = "NA", col7 = "NA", col8 = "NA"
  ) %>%
  dplyr::select(
    id, col1, col2, col3, TranscriptstableID,
    col4, VEP, col5, cdspos, protpos, aachange2, codonsub2
  )

dummy_ssb %<>% dplyr::mutate(aachange2 = as.character(aachange2))
dummy_ssb %<>% dplyr::mutate(codonsub2 = as.character(codonsub2))

df_ssb <- as.data.frame(dummy_ssb)

df_ssb <- df_ssb[!is.na(df_ssb$TranscriptstableID), ]

write.table(
  df_ssb,
  file = out_name, quote = FALSE, sep = "\t",
  col.names = FALSE, row.names = FALSE
)
