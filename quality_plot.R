# Quality Plots ################################################################
# We use DADA2 to obtain quality plots which we will use to filter in a later
# section

## Load Libraries ==============================================================
# Load all R packages you may need if not coming directly from the previous
# step.
library(dada2)
library(digest)
library(tidyverse)
library(seqinr)
library(ShortRead)

## File Housekeeping ===========================================================

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.

args <- commandArgs(trailingOnly = TRUE)

# Save project name as an object
project_name <- basename(dirname(getwd()))
print("This project is named")
project_name
# Set a path to the directory with the cutadapt-trimmed reads.
path_to_trimmed <- "../data/working/trimmed_sequences"

trimmed_F <- sort(
  list.files(
    path_to_trimmed,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
trimmed_R <- sort(
  list.files(
    path_to_trimmed,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)

# Make a new vector of sample names from your trimmed reads.
sample_names_trimmed <- sapply(
  strsplit(basename(trimmed_F), "_S\\d{1,3}_"),
  `[`,
  1
)

## Get Read Counts of Raw Samples ==============================================
# Make a list of all the files in your "data/raw" folder.
reads_to_trim <- list.files("../data/raw")

# Separate files by read direction (R1,R2), and save each
reads_to_trim_F <- reads_to_trim[str_detect(reads_to_trim, "R1_001.fastq.gz")]

# Separate the elements of "reads_to_trim_F" by underscore, and save the first
# element as "sample_names".
sample_names_raw <- sapply(
  strsplit(basename(reads_to_trim_F), "_S\\d{1,3}_"),
  `[`,
  1
)

# Count the number of reads in each sample.
sequence_counts_raw <- sapply(
  paste0("../data/raw/", reads_to_trim_F),
  function(file) {
    fastq_data <- readFastq(file)
    length(fastq_data)
  }
)
# Name these counts with your sample names
names(sequence_counts_raw) <- sample_names_raw

## Get Read Counts of Trimmed Samples ==========================================
# Count the number of reads in each trimmed sample. Since cutadapt only
# keeps paired reads, we only need to count forward samples.
sequence_counts_trimmed <- sapply(trimmed_F, function(file) {
  fastq_data <- readFastq(file)
  length(fastq_data)
})
names(sequence_counts_trimmed) <- sample_names_trimmed

## Remove empty sample files ===================================================
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_noreads_F <- trimmed_F[sapply(trimmed_F, file_size) < 100]
trimmed_noreads_R <- trimmed_R[sapply(trimmed_R, file_size) < 100]

print(
  "Here are the sample files for which contain no reads after primer trimming:"
)
names(trimmed_noreads_F)

# Remove the empty files. invisible() does it without a logical (TRUE/FALSE)
# message
invisible(file.remove(trimmed_noreads_R))
invisible(file.remove(trimmed_noreads_F))

# Renew vector of sample names from your trimmed reads.
sample_names_trimmed <- sapply(
  strsplit(basename(trimmed_F), "_S\\d{1,3}_"),
  `[`,
  1
)
## Make Quality Plots ==========================================================

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
quality_plot_F <- plotQualityProfile(
  trimmed_F[1:length(sample_names_trimmed)],
  aggregate = TRUE
)
quality_plot_F_reduced <- quality_plot_F +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )
# Examine the reverse reads as you did the forward.
quality_plot_R <- plotQualityProfile(
  fnRs[1:length(sample_names_trimmed)],
  aggregate = TRUE
)
quality_plot_R_reduced <- quality_plot_R +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )
# Export quality plots.
ggsave(
  paste0("../data/results/", project_name, "_qualplotF.pdf"),
  plot = quality_plot_F,
  width = 9,
  height = 9
)
ggsave(
  paste0("../data/results/", project_name, "_qualplotR.pdf"),
  plot = quality_plot_R,
  width = 9,
  height = 9
)

# Save all the objects created to this point in this section
save.image(file = "../data/working/1_trim_qual.RData")
