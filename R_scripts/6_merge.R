# MERGE ########################################################################
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
    stop(paste0("No gene argument provided in job file 6_merge_", gene, ".job"))
}

# get gene name from argument
gene <- args[1]

# Load the RData from the denoising step
load(paste0("data/working/5_denoise_", gene, ".RData"))

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. mergePairs calls for the forward denoising
# result (denoised), then the forward filtered and truncated reads
# (filtered_reads), then the same for the reverse reads.

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.

# Merge sequences with a 12 bp minimum overflap and no mismatches allowed in the overlap region. 
# This is a pretty strict merge, but I have found that it is better to be strict here and lose 
# some reads than to be more lenient and have more spurious sequences in the final results.
cat("\nMerging forward and reverse reads for", gene, "\n")
merged_reads <- mergePairs(
  denoised$F,
  filtered_reads$F,
  denoised$R,
  filtered_reads$R,
  minOverlap = 12,
  maxMismatch = 0,
  verbose = TRUE
)
cat("\nMerging is complete for", gene, "\n")

## Create Sequence-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.

seqtab <- makeSequenceTable(merged_reads)

# This describes the dimensions of the gene-specific table just made
# First the number of samples
cat(
  "\nThis is the number of samples for your",
  gene,
  "Sequence-Table:",
  length(rownames(seqtab)),
  "\n"
)
# Then the number of ASVs
cat(
  "\nThis is the number of ASVs for your",
  gene,
  "Sequence-Table:",
  length(colnames(seqtab)),
  "\n"
)

# Save all the objects created to this point in this section
save.image(file = paste0("data/working/6_merge_", gene, ".RData"))

cat("6_merge.job is complete and reads have been merged for", gene, "\n")
