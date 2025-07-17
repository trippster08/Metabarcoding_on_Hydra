# Quality Plots ################################################################
# We use DADA2 to obtain quality plots which we will use to filter in a later
# section

## Load Libraries ==============================================================
# Load all R packages you may need. suppressMessages keeps it exporting updates
# to the log
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Save project name as an object
project_name <- basename(getwd())
cat("\nThis project is named ", project_name, ".\n")

# Get list of genes passed from job file
genes <- commandArgs(trailingOnly = TRUE)
# Determine number of genes to analyze
num_genes <- length(genes)
# Prints to log a list of genes that will be analyzed
cat("\nWe will be creating quality plots for:", genes, "\n")

# Set a path to the directory containing raw reads.
path_to_raw_reads <- "data/raw"
# Set path to working directory
path_to_working <- "data/working"
# Set path to the directory (or directories, depending upon the number of genes)
# of the trimmed sequences. This creates a list of paths, one path for each gene
path_to_trimmed <- setNames(
  lapply(genes, function(gene) paste0("data/working/trimmed_sequences/", gene)),
  genes
)
# Set a path to the results directorie(s), as a list of paths, one for each gene
path_to_results <- setNames(
  lapply(genes, function(gene) paste0("data/results/", gene)),
  genes
)

## Get Raw Read Counts =========================================================
# Make a list of all the files in your "data/raw" folder.
reads_to_trim <- list.files(path_to_raw_reads)

# Get the names of all the samples in "data/raw". The name is everything
# before the illumina barcode
sample_names_raw <- sapply(
  strsplit(
    basename(
      reads_to_trim[str_detect(reads_to_trim, "R1_001.fastq.gz")]
    ),
    "_S\\d{1,3}_"
  ),
  `[`,
  1
)


# Count the number of reads in each sample.
sequence_counts_raw <- sapply(
  paste(
    path_to_raw_reads,
    reads_to_trim[str_detect(reads_to_trim, "R1_001.fastq.gz")],
    sep = "/"
  ),
  function(file) {
    fastq_data <- readFastq(file)
    length(fastq_data)
  }
)

# Name these counts with your sample names
names(sequence_counts_raw) <- sample_names_raw
# Print to log the number of reads for each sample
cat("\nHere are the raw read counts for each sample:\n")
print(sequence_counts_raw)

## Trimmed Reads ===============================================================
# Make a list of all gene-specific trimmed reads (it's a list of 3 gene-specific
# vectors containing trimmed read names), with each gene-specific item
# containing a list of Forward(R1) and a list of Reverse(R2) reads.
trimmed_reads <- setNames(
  lapply(genes, function(gene) {
    list(
      F = sort(list.files(
        path_to_trimmed[[gene]],
        pattern = "_R1.fastq.gz",
        full.names = TRUE
      )),
      R = sort(list.files(
        path_to_trimmed[[gene]],
        pattern = "_R2.fastq.gz",
        full.names = TRUE
      ))
    )
  }),
  genes
)
# Make list of gene-specific sample names for trimmed reads, and populate list
# with names from gene-specific trimmed read files.
sample_names_trimmed <- setNames(
  lapply(genes, function(gene) {
    files <- trimmed_reads[[gene]]$F
    if (length(files) == 0 || is.null(files)) {
      return(character(0))
    }
    sapply(
      strsplit(basename(files), "_trimmed"),
      `[`,
      1
    )
  }),
  genes
)

## Remove empty sample files ===================================================
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.

# Create a list to store non-empty trimmed reads
actual_trimmed_reads <- setNames(vector("list", length(genes)), genes)
# Create a lis to store trimmed read counts
sequence_counts_trimmed <- setNames(vector("list", length(genes)), genes)

# For each gene, create a list of trimmed forward and trimmed reverse reads.
# These are ephemeral for each gene, and are not
for (gene in genes) {
  # Count reads in R1 files
  sequence_counts_trimmed[[gene]] <- sapply(
    trimmed_reads[[gene]]$F,
    function(file) {
      fq <- readFastq(file)
      length(fq)
    }
  )

  # Name the read counts using sample_names_trimmed and store them
  names(sequence_counts_trimmed[[gene]]) <- sample_names_trimmed[[gene]]

  # Keep only files with at least 1 read
  valid_indices <- which(sequence_counts_trimmed[[gene]] > 0)

  actual_trimmed_reads[[gene]] <- list(
    F = trimmed_reads[[gene]]$F[valid_indices],
    R = trimmed_reads[[gene]]$R[valid_indices]
  )
}

# Remove read files with zero reads for each gene
for (gene in genes) {
  original_samples <- sample_names_trimmed[[gene]]
  kept_files <- actual_trimmed_reads[[gene]]$F

  # Extract sample names from kept files
  kept_samples <- if (length(kept_files) > 0) {
    sapply(strsplit(basename(kept_files), "_trimmed"), `[`, 1)
  } else {
    character(0)
  }

  # Reset names for trimmed without removed samples
  sample_names_trimmed <- setNames(
    lapply(genes, function(gene) {
      files <- actual_trimmed_reads[[gene]]$F
      if (length(files) == 0 || is.null(files)) {
        return(character(0))
      }

      sapply(
        strsplit(basename(files), "_trimmed"),
        `[`,
        1
      )
    }),
    genes
  )

  # Find filtered-out samples
  removed_samples <- setdiff(original_samples, kept_samples)

  # Print results
  if (length(removed_samples) > 0) {
    cat(
      "\nAfter trimming,",
      length(removed_samples),
      "samples had zero reads for",
      gene,
      "and were removed:\n"
    )
    print(removed_samples)
  } else {
    cat("\nAll samples had zero reads for", gene, "\n")
  }
}

# Count the number of samples remaining for each gene, and print
for (gene in genes) {
  nsamples <- length(sample_names_trimmed[[gene]])
  cat("\nWe will analyze", nsamples, "samples for", gene, "\n")
}


### Make Quality Plots ---------------------------------------------------------

# This creates DADA2 quality plots for each gene. It aggregates all samples and
# creates a single plot. These are the same quality plots as Qiime2, but
# maybe slightly easier to interpret.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.

# Create gene-specific list for the quality plots that will be saved and printed
quality_plots_better <- setNames(vector("list", length(genes)), genes)

# For each gene, create a list of F and R for storing the forward and reverse
# plots
for (gene in genes) {
  quality_plots <- list(F = NULL, R = NULL)
  quality_plots_better[[gene]] <- list(F = NULL, R = NULL)
  # For each direction for each gene, create quality plots
  for (direction in c("F", "R")) {
    reads <- actual_trimmed_reads[[gene]][[direction]]
    sample_names <- sample_names_trimmed[[gene]]
    quality_plots <- plotQualityProfile(
      reads[1:length(sample_names)],
      aggregate = TRUE
    )
    plot_build <- ggplot_build(quality_plots)
    max_x <- plot_build$layout$panel_params[[1]]$x.range[2]
    # Make the quality plots easier to interpret by changing the x-axis scale,
    # creating vertical lines every 10 bp to better determine quality scores at
    # length, add name of gene to plot
    quality_plots_better[[gene]][[direction]] <- quality_plots +
      scale_x_continuous(
        limits = c(0, max_x),
        breaks = seq(0, max_x, 10)
      ) +
      geom_vline(
        xintercept = seq(0, max_x, 10),
        color = "blue",
        linewidth = 0.25
      ) +
      annotate(
        "text",
        x = 30,
        y = 2,
        label = paste0(gene, " ", direction),
        size = 6
      )

    # Save the plot as a PDF
    ggsave(
      filename = file.path(
        path_to_working,
        paste0(
          project_name,
          "_",
          gene,
          "_qualplot",
          direction,
          ".pdf"
        )
      ),
      plot = quality_plots_better[[gene]][[direction]],
      width = 9,
      height = 9
    )
  }
}

save.image(file = "data/working/2_qual.RData")

cat("\n2_quality.job is complete and quality plots have been created.\n")
