# DADA2 ########################################################################
# We use Dada2 to filter and trim reads, estimate error rates and use these
# estimates to denoise reads, merge paired reads, and remove chimeric sequences

## Load Libraries ==============================================================
# Load all R packages you may need if not coming directly from the previous
# step.
library(dada2)
library(digest)
library(tidyverse)
library(seqinr)
library(ShortRead)

## File Housekeeping ===========================================================

args <- commandArgs(trailingOnly = TRUE)
gene1 <- args[1]
gene2 <- args[2]

# Save project name as an object
project_name <- basename(getwd())

# Set a path to the directory with the cutadapt-trimmed reads.
path_to_trimmed_gene1 <- paste0("../data/working/trimmed_reads/", gene1)
path_to_trimmed_gene2 <- paste0("../data/working/trimmed_reads/", gene2)

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

# Make a new vector of sample names from your trimmed reads.
sample_names_trimmed_gene1 <- sapply(
  strsplit(basename(trimmed_gene1_F), "_"),
  `[`,
  1
)

sample_names_trimmed_gene2 <- sapply(
  strsplit(basename(trimmed_gene2_F), "_"),
  `[`,
  1
)

## Get Read Counts of Trimmed Samples ==========================================
# Count the number of reads in each trimmed sample. Since cutadapt only
# keeps paired reads, we only need to count forward samples.
sequence_counts_trimmed_gene1 <- sapply(trimmed_gene1_F, function(file) {
  fastq_data_gene1 <- readFastq(file)
  length(fastq_data_gene1)
})

sequence_counts_trimmed_gene2 <- sapply(trimmed_gene2_F, function(file) {
  fastq_data_gene2 <- readFastq(file)
  length(fastq_data_gene2)
})
names(sequence_counts_trimmed_gene1) <- sample_names_trimmed_gene1
names(sequence_counts_trimmed_gene2) <- sample_names_trimmed_gene2


## Remove empty sample files ===================================================
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_noreads_gene1_F <- trimmed_gene1_F[
  sapply(trimmed_gene1_F, file.size) < 100
]
file.remove(trimmed_noreads_gene1_F)
trimmed_noreads_gene1_R <- trimmed_gene1_R[
  sapply(trimmed_gene1_R, file.size) < 100
]
file.remove(trimmed_noreads_gene1_R)

print(paste(
  "Here are the samples files for",
  gene1,
  "which contain no reads after primer trimming ",
  sep = " "
))
trimmed_noreads_gene1_F
# This saves the gene2 R1 and R2 fastq sample files only if both the R1 and R2
# sample files have reads.
trimmed_noreads_gene2_F <- trimmed_gene2_F[
  sapply(trimmed_gene2_F, file.size) < 100
]
file.remove(trimmed_noreads_gene2_F)
trimmed_noreads_gene2_R <- trimmed_gene2_R[
  sapply(trimmed_gene2_R, file.size) < 100
]
file.remove(trimmed_noreads_gene2_R)

print(paste(
  "Here are the samples files for",
  gene2,
  "which contain no reads after primer trimming ",
  sep = " "
))
trimmed_noreads_gene2_F

### Remove mistmatched reads----------------------------------------------------
# Check to see how many wrong-gene occurances there are for gene1 and save the
# names of the samples with these mismatches.
mismatches_gene1 <- sort(
  list.files(
    trimmed_gene1,
    pattern = gene2,
    full.names = TRUE
  )
)

# Move all the misidentified/empty files into a newly created "mismatches"
# directory. file.rename moves the files you want to move, and deletes them from
# their original directory.
file.rename(
  from = mismatches_gene1,
  to = paste0(
    "../data/working/trimmed_reads/mismatches/",
    basename(mismatches_gene1)
  )
)

# Repeat this process with your second gene. Make sure to reverse the path to
# your trimmed reads, and "pattern=" arguments

mismatches_gene2 <- sort(
  list.files(
    paste0("../data/working/trimmed_reads/", gene2),
    pattern = gene1,
    full.names = TRUE
  )
)

# Move all the misidentified/empty files into a newly created "misID_gene1"
# directory. file.rename moves the files you want to move, and deletes them from
# their original directory.
file.rename(
  mismatches_gene2,
  to = paste0(
    "../data/working/trimmed_reads/mismatches/",
    basename(mismatches_gene2)
  )
)

## Gene1 =======================================================================

nsamples_gene1 <- length(sample_names_gene1)
print(paste(
  "Here are the number of",
  gene1,
  "samples that will be analyzed with DADA2: ",
  sep = " "
))
nsamples_gene1

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
    "/qualplotF_",
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
    "/qualplotR_",
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

nsamples_gene2 <- length(sample_names_gene2)

print(paste(
  "Here are the number of",
  gene1,
  "samples that will be analyzed with DADA2: ",
  sep = " "
))
nsamples_gene2

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
    "/qualplotF_",
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
    "/qualplotR_",
    gene2,
    ".pdf"
  ),
  plot = quality_plot_gene2_R,
  width = 9,
  height = 9
)

save(
  gene1,
  gene2,
  path_to_trimmed_gene1,
  path_to_trimmed_gene2,
  trimmed_gene1_F,
  trimmed_gene1_R,
  trimmed_gene2_F,
  trimmed_gene2_R,
  sample_names_trimmed_gene1,
  sample_names_trimmed_gene2,
  sequence_counts_trimmed_gene1,
  sequence_counts_trimmed_gene2,
  trimmed_noreads_gene1_F,
  trimmed_noreads_gene1_R,
  trimmed_noreads_gene2_F,
  trimmed_noreads_gene2_R,
  mismatches_gene1,
  mismatches_gene2,
  nsamples_gene1,
  nsamples_gene2,
  quality_plot_gene1_F,
  quality_plot_gene1_F_reduced,
  quality_plot_gene1_R,
  quality_plot_gene1_R_reduced,
  quality_plot_gene2_F,
  quality_plot_gene2_F_reduced,
  quality_plot_gene2_R,
  quality_plot_gene2_R_reduced,
  file = "data/working/1_trim_qual.Rdata"
)
