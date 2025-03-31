# FILTER DENOISE MERGE #########################################################
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## Merge Paired Sequences ======================================================

# Load the RData from "quality_plot_multigene.R"
load("../data/working/5_denoise.RData")

# Here we merge the paired reads. merged calls for the forward denoising result
# (denoised_F), then the forward filtered and truncated reads (filtered_F),
# then the same for the reverse reads (denoised_R and filtered_R).

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
save.image(file = "../data/working/6_merge.RData")

print("Job 6_merge.job has finished")
