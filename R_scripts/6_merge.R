# MERGE ########################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================

# Load the RData from "quality_plot_multigene.R"
load("data/working/5_denoise.RData")

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

# First, make a list to hold the gene-specific merged reads.
merged_reads <- setNames(vector("list", length(genes)), genes)
# Loop through each gene, creating gene-specific merged sequences
for (gene in genes) {
  cat("\nMerging forward and reverse reads for", gene)
  merged_reads[[gene]] <- mergePairs(
    denoised[[gene]]$F,
    filtered_reads[[gene]]$F,
    denoised[[gene]]$R,
    filtered_reads[[gene]]$R,
    minOverlap = 12,
    maxMismatch = 0,
    verbose = TRUE
  )
  cat("\nMerging is complete for", gene)
}
## Create Sequence-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.

# First make a list to hold the gene-specific sequence-table
seqtab <- setNames(vector("list", length(genes)), genes)
# Loop through each gene, making the sequence-table
for (gene in genes) {
  seqtab[[gene]] <- makeSequenceTable(merged_reads[[gene]])

  # This describes the dimensions of the gene-specific table just made
  # First the number of samples
  cat(
    "\nThis is the number of samples for your",
    gene,
    "Sequence-Table:",
    length(rownames(seqtab[[gene]])),
    "\n"
  )
  # Then the number of ASVs
  cat(
    "\nThis is the number of ASVs for your",
    gene,
    "Sequence-Table:",
    length(colnames(seqtab[[gene]])),
    "\n"
  )
}
# Save all the objects created to this point in this section
save.image(file = "data/working/6_merge.RData")

cat("6_merge.job is complete and reads have been merged.\n")
