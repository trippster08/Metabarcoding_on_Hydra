# DENOISE ######################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================

# Load the RData from "quality_plot_multigene.R"
load("data/working/4_error.RData")

## Denoising ===================================================================

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtered_F), and "err =" which is the error file from
# learnErrors (effF).

# First make a list to store the gene-specific denoised sequences
denoised <- setNames(vector("list", length(genes)), genes)
# Loop through all genes, denoising reads, after adding read direction to the
# list.
for (gene in genes) {
  denoised[[gene]] <- list(F = NULL, R = NULL)
  for (direction in c("F", "R")) {
    denoised[[gene]][[direction]] <- dada(
      filtered_reads[[gene]][[direction]],
      err = errors[[gene]][[direction]],
      errorEstimationFunction = loessErrfun,
      selfConsist = FALSE,
      pool = FALSE,
      multithread = TRUE,
      verbose = TRUE
    )
  }
  cat("Denoising is complete for", gene)
}

# Save all the objects created to this point in this section
save.image(file = "data/working/5_denoise.RData")

cat("\nJob 5_denoise.job is finished and data has been denoised by DADA2")
