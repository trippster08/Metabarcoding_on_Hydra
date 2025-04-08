# MERGE ########################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
args <- commandArgs(trailingOnly = TRUE)
gene <- args[1]

# Load the RData from "quality_plot_multigene.R"
load(paste0("../data/working/", gene, "_5_denoise.RData"))

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. merged calls for the forward denoising result
# (dadaFs), then the forward filtered and truncated reads (filtFs), then the
# same for the reverse reads (dadaRs and filtRs).

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.
merged_reads <- mergePairs(
  denoised_F,
  filtered_F,
  denoised_R,
  filtered_R,
  minOverlap = 12,
  maxMismatch = 0,
  verbose = TRUE
)

## Create Sequence-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab <- makeSequenceTable(merged_reads)
# This describes the dimensions of the table just made

print(paste(
  "This is the number of samples for your Sequence-Table:",
  length(rownames(seqtab)),
  sep = " "
))

print(paste(
  "This is the number of ASVs for your Sequence-Table:",
  length(colnames(seqtab)),
  sep = " "
))

# Save all the objects created to this point in this section
save.image(file = paste0("../data/working/", gene, "_6_merge.RData"))

print("Job 6_merge_multigene.job and this analysis has finished")
