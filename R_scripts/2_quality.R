# Quality Plots ################################################################
# We use DADA2 to obtain quality plots which we will use to filter in a later
# section

## Load Libraries ==============================================================
# Load all R packages you may need.
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
genes <- commandArgs(trailingOnly = TRUE)
num_genes <- length(genes)
gene_list <- setNames(as.list(genes), genes)
cat(
  "These are the genes we will be creating quality plots for: ",
  paste(genes, collapse = ", "),
  "\n"
)

# Save project name as an object
project_name <- basename(getwd())
print(paste0("This project is named ", project_name))

# Set a path to the directory with the raw reads and cutadapt-trimmed reads.
path_to_raw_reads <- "data/raw"
path_to_trimmed <- setNames(
  lapply(genes, function(gene) paste0("data/working/trimmed_sequences/", gene)),
  genes
)
path_to_results <- "data/results"
## Get Raw Read Counts =========================================================
# Make a list of all the files in your "data/raw" folder.
reads_to_trim <- list.files(path_to_raw_reads)

# Get the names of all the samples in "data/raw"
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

print(paste0("Here are the raw read counts for each sample:"))
sequence_counts_raw

## Trimmed Reads ===============================================================
# Create vectors for the trimmed reads, both forward (R1) and reverse (R2) and
# for gene1 and gene2.
# Loop through each gene name

# Make a list of all gene-specific trimmed reads (it's a list of 3 gene-specific
# vectors containing trimmed read names)

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

# Create a new list to store non-empty trimmed reads
actual_trimmed_reads <- setNames(vector("list", length(genes)), genes)
sequence_counts_trimmed <- setNames(vector("list", length(genes)), genes)


for (gene in genes) {
  F_reads <- trimmed_reads[[gene]]$F
  R_reads <- trimmed_reads[[gene]]$R
  sample_names <- sample_names_trimmed[[gene]]

  # Count reads in R1 files
  read_counts <- sapply(F_reads, function(file) {
    fq <- readFastq(file)
    length(fq)
  })

  # Name the read counts using sample_names_trimmed and store them
  names(read_counts) <- sample_names
  sequence_counts_trimmed[[gene]] <- read_counts

  # Keep only files with at least 1 read
  valid_indices <- which(read_counts > 0)

  actual_trimmed_reads[[gene]] <- list(
    F = F_reads[valid_indices],
    R = R_reads[valid_indices]
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
      "\nThese samples had zero reads for this gene and were removed:",
      gene,
      "\n"
    )
    print(removed_samples)
  } else {
    cat("\nNo samples had zero reads for this gene:", gene, "\n")
  }
}

#Count the number of samples remaining for each gene, and print
for (gene in genes) {
  nsamples <- length(sample_names_trimmed[[gene]])
  message <- paste("We will analyze", nsamples, "samples for", gene)
  print(message)
}


### Make Quality Plots ---------------------------------------------------------

# This visualizes the quality plots. If you want to look at quality plots for
# each individual sample, use "aggregate = FALSE", and include whichever sample
# number you want in the square brackets (to aggregate all samples, replace N
# with the number of samples, or with length(fnFs)). For example, "fnFs[1:2]"
# will result in two plots, one for the first sample and one for the second.
# "fnFs[1:17]" will result in 17 plots, one for each sample.  Using
# "aggregate = TRUE" will combine any samples called (for example, "fnFS[1:17]"
# aggregates sample 1 through 17) into a single plot. This results in the same
# the quality plots as Qiime2.

# For these plots, the green line is the mean quality score at that position,
# the orange lines are the quartiles (solid for median, dashed for 25% and 75%)
# and the red line represents the proportion of reads existing at that position.

quality_plots <- setNames(vector("list", length(genes)), genes)
quality_plots_better <- setNames(vector("list", length(genes)), genes)

for (gene in genes) {
  quality_plots[[gene]] <- list(F = NULL, R = NULL)
  quality_plots_better[[gene]] <- list(F = NULL, R = NULL)

  for (direction in c("F", "R")) {
    reads <- actual_trimmed_reads[[gene]][[direction]]
    sample_names <- sample_names_trimmed[[gene]]
    quality_plots <- plotQualityProfile(
      reads[1:length(sample_names)],
      aggregate = TRUE
    )
    plot_build <- ggplot_build(quality_plots)
    max_x <- plot_build$layout$panel_params[[1]]$x.range[2]

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
        path_to_results,
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

print("Job 2_quality.job and this analysis has finished")
