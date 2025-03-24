# FILTER DENOISE MERGE #########################################################
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## Denoising ===================================================================

# Load the RData from "quality_plot_multigene.R"
load("../data/working/4_error.RData")

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtered_F), and "err =" which is the error file from
# learnErrors (effF).
denoised_F <- dada(
  filtered_F,
  err = errors_F,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

denoised_R <- dada(
  filtered_R,
  err = errors_R,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# Save all the objects created to this point in this section
save.image(file = "../data/working/5_denoise.RData")
