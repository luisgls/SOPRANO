# Load necessary libraries
suppressMessages(library("reshape"))
suppressMessages(library("tidyr"))

# Input data sources
args <- commandArgs(trailingOnly = TRUE)
file.name <- args[1]
file.name2 <- args[2]
file.name3 <- args[3]
file.name4 <- args[4]

# Read input data
df <- read.csv(file.name, header = FALSE, sep = "\t")
colnames(df) <- c("EnsembleID", "Total", "Class")

# Transform data from long to wide format
df2 <- spread(df, Class, Total, fill = 0)
df2 <- as.data.frame(df2)

# Read calculation for total sites
df.sites.extra <- read.csv(file.name2, header = FALSE, sep = "\t")
df.sites.intra <- read.csv(file.name3, header = FALSE, sep = "\t")

# Read intronic data
df.intron <- read.csv(file.name4, header = FALSE, sep = "\t")
colnames(df.intron) <- c("EnsemblID", "intronrate", "mutsintron", "intronlength")

# Merge data frames
df3 <- merge(df2, df.intron, by.x = "EnsembleID", by.y = "EnsemblID", all.x = TRUE)
df3[is.na(df3)] <- 0

# Calculate sums
na_epi <- sum(df3$extra_missense_variant)
ns_epi <- sum(df3$extra_synonymous_variant)
na_nonepi <- sum(df3$intra_missense_variant)
ns_nonepi <- sum(df3$intra_synonymous_variant)
ns_intron <- sum(df3$mutsintron)
NS_intron <- sum(df3$intronlength)
NS2_intron <- ns_intron / ((ns_nonepi / df.sites.intra[1, 2] + ns_epi / df.sites.extra[1, 2]) / 2)

# Calculate Ka/Ks ratios
ALL.extraKaKs <- (na_epi / df.sites.extra[1, 1]) / (ns_epi / df.sites.extra[1, 2])
ALL.intraKaKs <- (na_nonepi / df.sites.intra[1, 1]) / (ns_nonepi / df.sites.intra[1, 2])
ALL.intra_intronKaKs <- (na_nonepi / df.sites.intra[1, 1]) / ((ns_nonepi + ns_intron) / (df.sites.intra[1, 2] + NS2_intron))

# Function to calculate confidence intervals and p-values
calculate_ci_pval <- function(m, s, M, S) {
  p1 <- m / (M + 1)
  p2 <- s / (S + 1)
  globaldnds <- p1 / p2
  N1 <- M
  N2 <- S
  SE <- sqrt((1 - p1) / (N1 * p1) + (1 - p2) / (N2 * p2))
  finalLowCI <- globaldnds * exp(-1.96 * SE)
  finalHighCI <- globaldnds * exp(1.96 * SE)
  N <- m + s
  return(list(finalLowCI = finalLowCI, finalHighCI = finalHighCI, N = N, globaldnds = globaldnds))
}

# Calculate confidence intervals and p-values for epitope and non-epitope regions
ci_epi <- calculate_ci_pval(na_epi, ns_epi, df.sites.extra[1, 1], df.sites.extra[1, 2])
ci_nonepi <- calculate_ci_pval(na_nonepi, ns_nonepi, df.sites.intra[1, 1], df.sites.intra[1, 2])

# Estimate P-value from confidence interval of the non-target region
high <- ci_nonepi$finalHighCI
low <- ci_nonepi$finalLowCI
val <- ALL.extraKaKs - ALL.intraKaKs
SE <- (high - low) / (2 * 1.96)
EST <- val
z <- EST / SE
PVAL1 <- exp(-0.717 * z - 0.416 * z^2)
PVAL2 <- exp(-0.717 * -z - 0.416 * -z^2)

PVAL <- if (PVAL2 > 0 & PVAL2 <= 1) PVAL2 else PVAL1
PVAL <- if (PVAL < 0.0001) 0.0001 else PVAL

# Print the results
cat(paste(paste("coverage", "ON_dnds", "ON_lowci", "ON_highci", "ON_muts",
                "OFF_dnds", "OFF_lowci", "OFF_highci", "OFF_muts",
                "Pval",
                "ON_na", "ON_NA", "ON_ns", "ON_NS",
                "OFF_na", "OFF_NA", "OFF_ns", "OFF_NS", "\t"), "\n"))

cat(paste(paste("ExonicOnly", ALL.extraKaKs, ci_epi$finalLowCI, ci_epi$finalHighCI, ci_epi$N,
                ALL.intraKaKs, ci_nonepi$finalLowCI, ci_nonepi$finalHighCI, ci_nonepi$N,
                PVAL,
                na_epi, df.sites.extra[1, 1], ns_epi, df.sites.extra[1, 2],
                na_nonepi, df.sites.intra[1, 1], ns_nonepi, df.sites.intra[1, 2], "\t"), "\n"))

# Calculate confidence intervals and p-values for intronic regions
ci_intron <- calculate_ci_pval(na_nonepi, ns_nonepi + ns_intron, df.sites.intra[1, 1], df.sites.intra[1, 2] + NS2_intron)

# Estimate P-value from confidence interval of the non-target region with intronic data
high <- ci_intron$finalHighCI
low <- ci_intron$finalLowCI
val <- ALL.extraKaKs - ALL.intra_intronKaKs
SE <- (high - low) / (2 * 1.96)
EST <- val
z <- EST / SE
PVAL1intron <- exp(-0.717 * z - 0.416 * z^2)
PVAL2intron <- exp(-0.717 * -z - 0.416 * -z^2)

PVALintron <- if (PVAL2intron > 0 & PVAL2intron <= 1) PVAL2intron else PVAL1intron
PVALintron <- if (PVALintron < 0.0001) 0.0001 else PVALintron

# Print results with intronic data
cat(paste(paste("ExonicIntronic", ALL.extraKaKs, ci_epi$finalLowCI, ci_epi$finalHighCI, ci_epi$N,
                ALL.intra_intronKaKs, ci_intron$finalLowCI, ci_intron$finalHighCI, ci_intron$N,
                PVALintron,
                na_epi, df.sites.extra[1, 1], ns_epi, df.sites.extra[1, 2],
                na_nonepi, df.sites.intra[1, 1], (ns_nonepi + ns_intron), (df.sites.intra[1, 2] + NS2_intron),
                "\t"), "\n"))