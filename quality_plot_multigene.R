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

args <- commandArgs(trailingOnly = TRUE)
gene1 <- args[1]
gene2 <- args[2]

# Save project name as an object
project_name <- basename(dirname(getwd()))
print(paste0("This project is named ", project_name))

# Set a path to the directory with the raw reads and cutadapt-trimmed reads.
path_to_raw_reads <- "../data/raw"
path_to_trimmed_gene1 <- paste0("../data/working/trimmed_sequences/", gene1)
path_to_trimmed_gene2 <- paste0("../data/working/trimmed_sequences/", gene2)

## Get Raw Read Counts =========================================================
# Make a list of all the files in your "data/raw" folder.
reads_to_trim <- list.files(path_to_raw_reads)

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
  paste(path_to_raw_reads, reads_to_trim_F, sep = "/"),
  function(file) {
    fastq_data <- readFastq(file)
    length(fastq_data)
  }
)

# Name these counts with your sample names
names(sequence_counts_raw) <- sample_names_raw

print(paste0("Here are the read counts for each raw sample:"))
sequence_counts_raw

## Remove Mistmatched Trimmed Read Files =======================================
# Check to see how many wrong-gene occurances there are for gene1 and save the
# names of the samples with these mismatches.
mismatches_gene1 <- sort(
  list.files(
    path_to_trimmed_gene1,
    pattern = gene2,
    full.names = TRUE
  )
)

# Move all the misidentified/empty files into a newly created "mismatches"
# directory. file.rename moves the files you want to move, and deletes them from
# their original directory.
invisible(file.rename(
  from = mismatches_gene1,
  to = paste0(
    "../data/working/trimmed_sequences/mismatches/",
    basename(mismatches_gene1)
  )
))

# Repeat this process with your second gene. Make sure to reverse the path to
# your trimmed reads, and "pattern=" arguments

mismatches_gene2 <- sort(
  list.files(
    paste0("../data/working/trimmed_sequences/", gene2),
    pattern = gene1,
    full.names = TRUE
  )
)

# Move all the misidentified/empty files into a newly created "misID_gene1"
# directory. file.rename moves the files you want to move, and deletes them from
# their original directory.
invisible(file.rename(
  mismatches_gene2,
  to = paste0(
    "../data/working/trimmed_sequences/mismatches/",
    basename(mismatches_gene2)
  )
))

## Create vectors for the trimmed reads, both forward (R1) and reverse (R2) and
# for gene1 and gene2.
trimmed_gene1_F <- sort(
  list.files(
    path_to_trimmed_gene1,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
trimmed_gene1_R <- sort(
  list.files(
    path_to_trimmed_gene1,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)

trimmed_gene2_F <- sort(
  list.files(
    path_to_trimmed_gene2,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
trimmed_gene2_R <- sort(
  list.files(
    path_to_trimmed_gene2,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)

# Make a vector of sample names from your trimmed reads.
sample_names_trimmed_gene1 <- sapply(
  strsplit(basename(trimmed_gene1_F), "_S\\d{1,3}_"),
  `[`,
  1
)

sample_names_trimmed_gene2 <- sapply(
  strsplit(basename(trimmed_gene2_F), "_S\\d{1,3}_"),
  `[`,
  1
)
## Remove empty sample files ===================================================
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_noreads_gene1_F <- trimmed_gene1_F[
  sapply(trimmed_gene1_F, file.size) < 100
]
trimmed_noreads_gene1_R <- trimmed_gene1_R[
  sapply(trimmed_gene1_R, file.size) < 100
]

sample_names_trimmed_noreads_gene1 <- sapply(
  strsplit(basename(trimmed_noreads_gene1_F), "_S\\d{1,3}_"),
  `[`,
  1
)
names(trimmed_noreads_gene1_F) <- sample_names_trimmed_noreads_gene1

print(paste(
  "Here are the samples files for",
  gene1,
  "which contain no reads after primer trimming:",
  sep = " "
))
names(trimmed_noreads_gene1_F)

invisible(file.remove(trimmed_noreads_gene1_F))
invisible(file.remove(trimmed_noreads_gene1_R))

# Remove zero read samples from trimmed_gene1_F amd trimmed_gene1_R vectors
trimmed_exists_gene1_F <- trimmed_gene1_F[
  !trimmed_gene1_F %in% trimmed_noreads_gene1_F
]
trimmed_exists_gene1_R <- trimmed_gene1_R[
  !trimmed_gene1_R %in% trimmed_noreads_gene1_R
]
trimmed_gene1_F <- trimmed_exists_gene1_F
trimmed_gene1_R <- trimmed_exists_gene1_R

# Make a new vector of sample names from your trimmed reads.
sample_names_trimmed_gene1 <- sapply(
  strsplit(basename(trimmed_gene1_F), "_S\\d{1,3}_"),
  `[`,
  1
)

# This saves the gene2 R1 and R2 fastq sample files only if both the R1 and R2
# sample files have reads.
trimmed_noreads_gene2_F <- trimmed_gene2_F[
  sapply(trimmed_gene2_F, file.size) < 100
]
trimmed_noreads_gene2_R <- trimmed_gene2_R[
  sapply(trimmed_gene2_R, file.size) < 100
]

sample_names_trimmed_noreads_gene2 <- sapply(
  strsplit(basename(trimmed_noreads_gene2_F), "_S\\d{1,3}_"),
  `[`,
  1
)
names(trimmed_noreads_gene2_F) <- sample_names_trimmed_noreads_gene2
print(paste(
  "Here are the samples files for",
  gene2,
  "which contain no reads after primer trimming ",
  sep = " "
))
names(trimmed_noreads_gene2_F)

invisible(file.remove(trimmed_noreads_gene2_F))
invisible(file.remove(trimmed_noreads_gene2_R))

# Remove zero read samples from trimmed_gene2_F amd trimmed_gene2_R vectors
trimmed_exists_gene2_F <- trimmed_gene2_F[
  !trimmed_gene2_F %in% trimmed_noreads_gene2_F
]
trimmed_exists_gene2_R <- trimmed_gene2_R[
  !trimmed_gene2_R %in% trimmed_noreads_gene2_R
]
trimmed_gene2_F <- trimmed_exists_gene2_F
trimmed_gene2_R <- trimmed_exists_gene2_R

# Make a new vector of sample names from your trimmed reads.
sample_names_trimmed_gene2 <- sapply(
  strsplit(basename(trimmed_gene2_F), "_S\\d{1,3}_"),
  `[`,
  1
)

## Gene1 =======================================================================

nsamples_gene1 <- length(sample_names_trimmed_gene1)
print(paste(
  "We will analyze",
  nsamples_gene1,
  "samples for",
  gene1,
  sep = " "
))


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

quality_plot_gene1_F <- plotQualityProfile(
  trimmed_gene1_F[1:length(sample_names_trimmed_gene1)],
  aggregate = TRUE
)
quality_plot_gene1_F_reduced <- quality_plot_gene1_F +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

# Examine the reverse reads as you did the forward.
quality_plot_gene1_R <- plotQualityProfile(
  trimmed_gene1_R[1:length(sample_names_trimmed_gene1)],
  aggregate = TRUE
)
quality_plot_gene1_R_reduced <- quality_plot_gene1_R +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

### Export Quality Plots -------------------------------------------------------

ggsave(
  paste0(
    "../data/results/",
    project_name,
    "_",
    gene1,
    "_qualplotF_",
    gene1,
    ".pdf"
  ),
  plot = quality_plot_gene1_F,
  width = 9,
  height = 9
)

ggsave(
  paste0(
    "../data/results/",
    project_name,
    "_",
    gene1,
    "_qualplotR_",
    gene1,
    ".pdf"
  ),
  plot = quality_plot_gene1_R,
  width = 9,
  height = 9
)

## Gene2 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).

nsamples_gene2 <- length(sample_names_trimmed_gene2)

print(paste(
  "We will analyze",
  nsamples_gene2,
  "samples for",
  gene2,
  sep = " "
))

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
quality_plot_gene2_F <- plotQualityProfile(
  trimmed_gene2_F[1:length(sample_names_trimmed_gene2)],
  aggregate = TRUE
)
quality_plot_gene2_F_reduced <- quality_plot_gene2_F +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

# Examine the reverse reads as you did the forward.
quality_plot_gene2_R <- plotQualityProfile(
  trimmed_gene2_R[1:length(sample_names_trimmed_gene2)],
  aggregate = TRUE
)
quality_plot_gene2_R_reduced <- quality_plot_gene2_R +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

### Export Quality Plots -------------------------------------------------------

ggsave(
  paste0(
    "../data/results/",
    project_name,
    "_",
    gene2,
    "_qualplotF_",
    gene2,
    ".pdf"
  ),
  plot = quality_plot_gene2_F,
  width = 9,
  height = 9
)

ggsave(
  paste0(
    "../data/results/",
    project_name,
    "_",
    gene2,
    "_qualplotR_",
    gene2,
    ".pdf"
  ),
  plot = quality_plot_gene2_R,
  width = 9,
  height = 9
)


# Here we replace "gene1" and "gene2" in each object with the actual name of
# gene1 and gende2, so each object will be associated with the correct data
# in later steps
all_objects <- ls()
filtered_gene1 <- grep("gene1", all_objects, value = TRUE)
filtered_gene2 <- grep("gene2", all_objects, value = TRUE)
objects_to_remove <- c()

for (obj in filtered_gene1) {
  new_name <- sub("gene1", gene1, obj)
  assign(new_name, get(obj))
  objects_to_remove <- c(objects_to_remove, obj)
}

for (obj in filtered_gene2) {
  new_name <- sub("gene2", gene2, obj)
  assign(new_name, get(obj))
  objects_to_remove <- c(objects_to_remove, obj)
}
rm(list = objects_to_remove)

save.image(file = "../data/working/2_qual.RData")
