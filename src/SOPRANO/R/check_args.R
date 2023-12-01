check_source_path <- function(source_path) {
  if (is.null(source_path)) {
    stop("VCF source(s) path not defined. Flag -s | --sources required")
  }
  if (!file.exists(source_path)) {
    stop(paste("Source(s) path does not exist:", source_path))
  }
}

check_output_path <- function(output_path) {
  if (!is.null(output_path)) {
    if (file.exists(output_path)) {
      stop(paste("This process would overwrite", output_path))
    }
  }
}

check_translator_dir <- function(translator_dir) {
  if (is.null(translator_dir)) {
    stop("Translator directory not defined. Flag -t | --t <dir path> required")
  }
  if (!dir.exists(translator_dir)) {
    stop(paste("Translator directory does not exits:", translator_dir))
  }
}

check_assembly_type <- function(assembly_type) {
  if (!assembly_type %in% c("GRCh38", "GRCh37")) {
    stop(paste(
      "Currently only supporting GRCh37 and GRCh38 assemblies, not", assembly_type
    ))
  }
}

check_auxiliary_paths <- function(dir_path, ...) {
  file_names <- list(...)
  for (file_name in file_names) {
    p <- file.path(dir_path, file_name)
    if (!file.exists(p)) {
      stop(paste("Auxiliary file not found:", p))
    }
  }
}
