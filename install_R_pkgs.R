package_list <- c(
  "devtools", "vcfR", "readr", "stringr", "tidyr", "optparse", "poilog",
  "plyr", "dplyr", "tibble", "magrittr", "MASS"
)

get_installed <- function() {
  return(installed.packages()[, "Package"])
}

for (pkg in package_list) {
  if (!pkg %in% get_installed()) {
    install.packages(
        pkg,
        character.only = TRUE,
        repos = "https://cloud.r-project.org/"
    )
  }
  library(pkg, character.only = TRUE)
}

gh_package_list <- c("im3sanger/dndscv")

for (gh_pkg in gh_package_list) {
  pkg <- str_split(gh_pkg, pattern = "/")[[1]][2]
  if (!pkg %in% get_installed()) {
    install_github(gh_pkg, character.only = TRUE, dependencies = FALSE)
  }
  library(pkg, character.only = TRUE)
}
