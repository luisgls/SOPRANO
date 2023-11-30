select_vcf_files <- function(source_path) {
  if (dir.exists(source_path)) {
    vcf_files <- Sys.glob(file.path(source_path, "*.vcf.gz"))
  }
  else {
    vcf_files <- c(source_path)
  }
}

extract_biallelic_indels <- function(vcf_path) {
  vcf_data <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  vcf_data <- vcfR::extract.indels(vcf_data)
  vcf_data <- vcf_data[vcfR::is.biallelic(vcf_data),]
  return(vcf_data)
}

extract_gt_tidy <- function(vcf) {
  vcfR::extract_gt_tidy(vcf)
}


fix_vcf_fields <- function(vcf) {
  vcf@fix %>% tibble::as_tibble()
}

merge_gt_tidy_with_fixed_fields <- function(gt_tidied_data, fixed_fields) {
  mapply(c, fixed_fields, gt_tidied_data, SIMPLIFY = FALSE)
}

reduce_tibble_id_to_chr_pos_ref_alt <- function(merged_data) {
  reduced <- tibble::as_tibble(cbind(merged_data$CHROM, merged_data$POS, merged_data$REF, merged_data$ALT))
  reduced <- dplyr::bind_rows(reduced, .id = "name")
  names(reduced) <- c("sampleID", "chr", "position", "ref_allele", "alt_allele")
  return(reduced)
}

correct_chr_patterns <- function(data) {
  data$chr <- gsub("chr", "", as.vector(data$chr))
}

add_transcript_and_variant_info <- function(dndscv_result) {
  annotation <- dplyr::as_tibble(dndscv_result$annotmuts)
  annotation <- annotation %>%
    dplyr::left_join(., transcriptlist, by = c("pid" = "ProteinstableID")) %>%
    dplyr::left_join(., variantlist, by = c("impact" = "REF")) %>%
    dplyr::select(
      chr, pos, ref, mut, TranscriptstableID, VEP, aachange, ntchange, codonsub
    )
  return(annotation)
}

update_annotation_formatting <- function(annotated_data) {
  annotated_data$codonsub2 <- gsub(">", "/", annotated_data$codonsub)
  annotated_data$aachange2 <- gsub("[0-9]+", "/", annotated_data$aachange)
  number <- "[0-9]+"
  annotated_data <- annotated_data %>%
    tidyr::unite(idtmp, c("chr", "pos"), sep = "_") %>%
    tidyr::unite(change, c("ref", "mut"), sep = "/") %>%
    tidyr::unite(id, c("idtmp", "change"), sep = "_") %>%
    dplyr::mutate(
      protpos = stringr::str_extract(annotated_data$aachange, number),
      cdspos = stringr::str_extract(annotated_data$ntchange, number)
    )
  return(annotated_data)
}

prepare_ssb_input_file <- function(processed_data, output_path) {
  ssb_daata <- processed_data %>%
    dplyr::mutate(
      col1 = "NA", col2 = "NA", col3 = "NA", col4 = "Transcript", col5 = "NA",
      col6 = "NA", col7 = "NA", col8 = "NA"
    ) %>%
    dplyr::select(
      id, col1, col2, col3, TranscriptstableID,
      col4, VEP, col5, cdspos, protpos, aachange2, codonsub2
    )

  ssb_daata %<>% dplyr::mutate(aachange2 = as.character(aachange2))
  ssb_daata %<>% dplyr::mutate(codonsub2 = as.character(codonsub2))

  ssb_daata <- as.data.frame(ssb_daata)
  ssb_daata <- ssb_daata[!is.na(df_ssb$TranscriptstableID),]

  write.table(
    ssb_daata,
    file = output_path, quote = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE
  )
}