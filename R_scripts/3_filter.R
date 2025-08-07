# FILTER #######################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Load the RData from "quality_plot_multigene.R"
load("data/working/2_qual.RData")
cat("\nThis project is named", project_name, "\n\n")
# Get arguments from job file. This should include the number of genes, gene
# names, and the R1 and R2 truncation values for each gene
args <- commandArgs(trailingOnly = TRUE)
# First argument is gene number
gene_num <- as.numeric(args[1])
# Save vector of genes (determined by gene number)
genes <- args[2:(gene_num + 1)]
# Save vector of truncation values (also determined by gene number)
truncation_values <- as.integer(args[(gene_num + 2):length(args)])
# Create a list of R1 and R2 truncation values for each gene
truncation_list <- setNames(
  lapply(seq_len(gene_num), function(i) {
    start <- (i - 1) * 2 + 1
    end <- start + 1
    truncation_values[start:end]
  }),
  genes
)

# This creates file paths for the reads that will be quality filtered with dada2
# in the next step and a list for storing filtered reads.
filtered_reads <- setNames(vector("list", length(genes)), genes)
path_to_filtered <- setNames(
  lapply(genes, function(gene) {
    paste0("data/working/filtered_sequences/", gene)
  }),
  genes
)
# For each gene, create a list for forward and reverse reads in filtered_reads,
# and create a list of filtered sample names from trimmed sample
# names, and add these names to the filtered_reads list,
for (gene in genes) {
  sample_names_filtered <- sample_names_trimmed[[gene]]
  filtered_reads[[gene]] <- list(F = NULL, R = NULL)
  for (direction in c("F", "R")) {
    filtered_reads[[gene]][[direction]] <- setNames(
      file.path(
        path_to_filtered[[gene]],
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
for (gene in genes) {
  filterAndTrim(
    actual_trimmed_reads[[gene]]$F,
    filtered_reads[[gene]]$F,
    actual_trimmed_reads[[gene]]$R,
    filtered_reads[[gene]]$R,
    truncLen = c(truncation_list[[gene]]),
    maxN = 0,
    maxEE = c(2, 2),
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE,
    verbose = TRUE
  )
}

# Create a new sample_names_filtered, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
for (gene in genes) {
  sample_names_filtered <- setNames(
    lapply(genes, function(gene) {
      files <- list.files(path_to_filtered[[gene]], pattern = "filt_F")
      sapply(strsplit(basename(files), "_filt_F"), `[`, 1)
    }),
    genes
  )
}
# Update filtered_reads, since some samples no longer have
# reads after filtering, and therefore no longer exist in the directory
for (gene in genes) {
  for (direction in c("F", "R")) {
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
}
# Create a list to contain the gene-specific filtered read counts
sequence_counts_filtered <- setNames(vector("list", length(genes)), genes)

# For each gene, count the number of forward reads, and add names to the read
# counts
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
  cat("\nHere are the samples removed after filtering for", gene, "\n")
  print(removed)
}

# Save all the objects created to this point in this section
save.image(file = "data/working/3_filter.RData")

cat("\n3_filter.job is finished and reads have been filtered by DADA2.\n")
