# FILTER #######################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Make objects to fill from job script
args <- commandArgs(trailingOnly = TRUE)
gene_num <- as.numeric(args[1])
genes <- args[2:(gene_num + 1)]
truncation_values <- as.integer(args[(gene_num + 2):length(args)])
truncation_list <- setNames(
  lapply(seq_len(gene_num), function(i) {
    start <- (i - 1) * 2 + 1
    end <- start + 1
    truncation_values[start:end]
  }),
  genes
)

# Load the RData from "quality_plot_multigene.R"
load("data/working/2_qual.RData")

# This creates file paths for the reads that will be quality filtered with dada2
# in the next step.
filtered_reads <- setNames(vector("list", length(genes)), genes)
path_to_filtered <- setNames(
  lapply(genes, function(gene) {
    paste0("data/working/filtered_sequences/", gene)
  }),
  genes
)
for (gene in genes) {
  sample_names_filtered <- sample_names_trimmed[[gene]]
  filtered_reads[[gene]] <- list(
    F = setNames(
      file.path(
        path_to_filtered,
        paste0(sample_names_filtered, "_filt_F.fastq.gz")
      ),
      sample_names_filtered
    ),
    R = setNames(
      file.path(
        path_to_filtered,
        paste0(sample_names_filtered, "_filt_R.fastq.gz")
      ),
      sample_names_filtered
    )
  )
}

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(i,j)" is how you tell Dada2 where to truncate all forward (i) and
# reverse (j) reads. Using "0" means reads will not be truncated.
# maxEE sets how many expected errors are allowed before a read is filtered out.

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 6,6) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.

# This step will create any new directories necessary (whatever path
# filtered_reads uses).

for (gene in genes) {
  filterAndTrim(
    actual_trimmed_reads[[gene]]$F,
    filtered_reads[[gene]]$F,
    actual_trimmed_reads[[gene]]$R,
    filtered_reads[[gene]]$R,
    truncLen = c(truncation_list[[gene]]),
    maxN = 0,
    maxEE = c(4, 4),
    rm.phix = TRUE,
    truncQ = 2,
    compress = TRUE,
    multithread = TRUE,
    verbose = TRUE
  )
}

# Save the reduced sample_names_filtered, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
sample_names_filtered <- setNames(
  lapply(genes, function(gene) {
    files <- list.files(path_to_filtered[[gene]], pattern = "filt_F")
    if (length(files) == 0 || is.null(files)) {
      return(character(0))
    }
    sapply(strsplit(basename(files), "_filt_F"), `[`, 1)
  }),
  genes
)

# Update filtered_reads, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
for (gene in genes) {
  F_filtered <- list.files(
    path_to_filtered[[gene]],
    pattern = "_filt_F",
    full.names = TRUE
  )
  R_filtered <- list.files(
    path_to_filtered[[gene]],
    pattern = "_filt_R",
    full.names = TRUE
  )

  filtered_reads[[gene]] <- list(
    F = setNames(
      list.files(
        path_to_filtered[[gene]],
        pattern = "filt_F",
        full.names = TRUE
      ),
      sapply(strsplit(basename(F_filtered), "_filt_F"), `[`, 1)
    ),
    R = setNames(
      list.files(
        path_to_filtered[[gene]],
        pattern = "filt_R",
        full.names = TRUE
      ),
      sapply(strsplit(basename(R_filtered), "_filt_R"), `[`, 1)
    )
  )
}

sequence_counts_filtered <- setNames(vector("list", length(genes)), genes)

for (gene in genes) {
  sequence_counts_filtered[[gene]] <- sapply(
    filtered_reads[[gene]]$F,
    function(file) {
      fq <- readFastq(file)
      length(fq)
    }
  )
  sample_names <- sapply(
    strsplit(basename(filtered_reads[[gene]]$F), "_filt_F"),
    `[`,
    1
  )
  names(sequence_counts_filtered[[gene]]) <- sample_names
}

# Report removed samples
for (gene in genes) {
  removed <- setdiff(
    sample_names_trimmed[[gene]],
    sample_names_filtered[[gene]]
  )
  cat("\nHere are the samples removed after filtering for", gene, ":\n")
  print(removed)
}

# Save all the objects created to this point in this section
save.image(file = "data/working/3_filter.RData")

print("Job 3_filter.job and this analysis has finished")
