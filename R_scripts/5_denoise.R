# DENOISE ######################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Get argument containing the gene name from job file.
args <- commandArgs(trailingOnly = TRUE)

# Check to make sure there is an argument for the gene name
if (length(args) < 1) {
    stop(paste0("No gene argument provided in job file 5_denoise_", gene, ".job."))
}

# get gene name from argument
gene <- args[1]
# Load the RData from "quality_plot_multigene.R"
load(paste0("data/working/4_error_", gene, ".RData"))

## Denoising ===================================================================

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtered_F), and "err =" which is the error file from
# learnErrors (effF).

cat("\nDenoising reads for", gene, "\n")
# Make a list for forward and reverse denoised sequences
denoised <- list(F = NULL, R = NULL)
# denoise forward and reverse reads
for (direction in c("F", "R")) {
  denoised[[direction]] <- dada(
    filtered_reads[[direction]],
    err = errors[[direction]],
    errorEstimationFunction = loessErrfun,
    selfConsist = FALSE,
    pool = FALSE,
    multithread = TRUE,
    verbose = TRUE
  )
}
cat("\nDenoising is complete for", gene, "\n")


# Save all the objects created to this point in this section
save.image(file = paste0("data/working/5_denoise_", gene, ".RData"))

cat("\n5_denoise.job for", gene, "is finished and data has been denoised by DADA2.\n")
