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
merged_reads <- setNames(vector("list", length(genes)), genes)
for (gene in genes) {
  merged_reads[[gene]] <- mergePairs(
    denoised[[gene]]$F,
    filtered_reads[[gene]]$F,
    denoised[[gene]]$R,
    filtered_reads[[gene]]$R,
    minOverlap = 12,
    maxMismatch = 0,
    verbose = TRUE
  )
}
## Create Sequence-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab <- setNames(vector("list", length(genes)), genes)
for (gene in genes) {
  seqtab[[gene]] <- makeSequenceTable(merged_reads[[gene]])

  # This describes the dimensions of the table just made

  print(paste(
    "This is the number of samples for your",
    gene,
    "Sequence-Table:",
    length(rownames(seqtab[[gene]])),
    sep = " "
  ))

  print(paste(
    "This is the number of ASVs for your",
    gene,
    "Sequence-Table:",
    length(colnames(seqtab[[gene]])),
    sep = " "
  ))
}
# Save all the objects created to this point in this section
save.image(file = "data/working/6_merge.RData")

print("Job 6_merge.job and this analysis has finished")
