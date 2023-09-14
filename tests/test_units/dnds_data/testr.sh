#!/bin/bash

script="/mnt/c/Users/kmarzouk/software/SOPRANO/src/SOPRANO/scripts/calculateKaKsEpiCorrected_CI_intron_V3.R"
prefix="TCGA-05-4396"

Rscript $script "$prefix.data.epitopes" "$prefix.epitopes.nans" "$prefix.intra_epitopes.nans" "$prefix.intron.rate"
