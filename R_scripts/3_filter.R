# FILTER #######################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Get arguments from job file. This should include the gene name
# and the R1 and R2 truncation values
args <- commandArgs(trailingOnly = TRUE)

# Check to make sure there are three arguments
if (length(args) < 3) {
    stop(paste0("Not enough arguments provided in filter job file. We need three arguments: gene name, R1 truncation value, and R2 truncation value."))
}

# Load the RData from "quality_plot_multigene.R"
load(paste0("data/working/2_qual.RData"))
cat("\nThis project is named", project_name, "\n")

# First argument is gene number
gene <- args[1]

# Save 2nd and 3rd arguments as r1 and r2 truncation values
r1_trunc <- as.numeric(args[2])
r2_trunc <- as.numeric(args[3])

cat("\nFiltering reads for gene:", gene,
    "with truncLen =", r1_trunc, r2_trunc, "\n")

# If path_to_results was stored as a list by gene, reduce to just this gene
if (exists("path_to_results") && is.list(path_to_results)) {
  path_to_results <- path_to_results[[gene]]
}


# This creates file paths for the reads that will be quality filtered with dada2
# in the next step and a list for storing filtered reads.
path_to_filtered <- paste0("data/working/filtered_sequences/", gene)

# Rename sample_names_trimmed to be gene-specific, because this is now a gene-specific script, and we will be saving gene-specific RData files, so we want to make sure the objects in those files are gene-specific as well.
sample_names_filtered <- sample_names_trimmed[[gene]]

# Create a list for forward and reverse reads in filtered_reads,
# and create a list of filtered sample names from trimmed sample
# names, and add these names to the filtered_reads list,
filtered_reads <- list(F = NULL, R = NULL)
for (direction in c("F", "R")) {
  filtered_reads[[direction]] <- setNames(
    file.path(
      path_to_filtered,
      paste0(
        sample_names_filtered,
        "_filt_",
        direction,
        ".fastq.gz"
      )
    ),
    sample_names_filtered
  )
}


# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 4,4) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.
filterAndTrim(
  actual_trimmed_reads[[gene]]$F,
  filtered_reads$F,
  actual_trimmed_reads[[gene]]$R,
  filtered_reads$R,
  truncLen = c(r1_trunc, r2_trunc),
  maxN = 0,
  maxEE = c(4, 4),
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE,
  verbose = TRUE
)


# Create a new sample_names_filtered, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
files_F <- list.files(path_to_filtered, pattern = "filt_F", full.names = TRUE)
samples_names_filtered <- sapply(strsplit(basename(files_F), "_filt_F"), `[`, 1)

# Update filtered_reads, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
F_filtered <- list.files(
  path_to_filtered,
  pattern = "_filt_F",
  full.names = TRUE
)
R_filtered <- list.files(
  path_to_filtered,
  pattern = "_filt_R",
  full.names = TRUE
)

filtered_reads <- list(
  F = setNames(
    list.files(
      path_to_filtered,
      pattern = "filt_F",
      full.names = TRUE
    ),
    sapply(strsplit(basename(F_filtered), "_filt_F"), `[`, 1)
  ),
  R = setNames(
    list.files(
      path_to_filtered,
      pattern = "filt_R",
      full.names = TRUE
    ),
    sapply(strsplit(basename(R_filtered), "_filt_R"), `[`, 1)
  )
)

# Count the number of forward reads, and add names to the read
# counts
sequence_counts_filtered <- sapply(
  filtered_reads$F,
  function(file) {
    fq <- readFastq(file)
    length(fq)
  }
)
sample_names <- sapply(
  strsplit(basename(filtered_reads$F), "_filt_F"),
  `[`,
  1
)
names(sequence_counts_filtered) <- sample_names


# Report removed samples
removed <- setdiff(
  sample_names_trimmed[[gene]],
  sample_names_filtered
)
cat("\nHere are the", gene, "samples removed after filtering\n")
print(removed)


# Save all the objects created to this point in this section
save.image(file = paste0("data/working/3_filter_", gene, ".RData"))

cat("\n3_filter.job for ", gene, " is finished and reads have been filtered by DADA2.\n")
