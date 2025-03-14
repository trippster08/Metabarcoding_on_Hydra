# DADA2 ########################################################################
# We use Dada2 to filter and trim reads, estimate error rates and use these
# estimates to denoise reads, merge paired reads, and remove chimeric sequences

## Load Libraries ==============================================================
# Load all R packages you may need if not coming directly from the previous
# step.
library(dada2)
library(tidyverse)
library(seqinr)

## File Housekeeping ===========================================================

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.

# Set a path to the directory with the cutadapt-trimmed reads.

args <- commandArgs(trailingOnly = TRUE)
gene1 <- args[1]
gene2 <- args[2]
trimmed_gene1 <- paste0("../data/working/trimmed_reads/", gene1)
trimmed_gene1_R1 <- sort(
  list.files(
    trimmed_gene1,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
trimmed_gene1_R2 <- sort(
  list.files(
    trimmed_gene1,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)
trimmed_gene2 <- paste0("../data/working/trimmed_reads/", gene2)
trimmed_gene2_R1 <- sort(
  list.files(
    trimmed_gene2,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
trimmed_gene2_R2 <- sort(
  list.files(
    trimmed_gene2,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)

# We will now run the rest of this section in multiple times, once for each
# gene present in the Illumina run. Again, replace each "gene1", "gene2",
# "gene3", etc with your specific gene name.

## Remove Empty or Misidentified Reads =========================================

# When cutadapt searches for primers, it removes any reads for which it cannot
# find primers to remove. Some samples (especially blanks) often end up with no
# reads left in the file, but cutadapt still creates a trimmed R1 and R2 file
# for that sample. DADA2 chokes on empty files, so we need to remove these.

# Also, when cutadapt moves each read to it's primer-specific directory, it
# sometimes misidentifies reads and places them into the wrong directory, and it
# creates files for the incorrect gene (e.g. creating a file for gene1 reads in
# the gene2 directory) but does not put any reads into that file (so it ends up
# being an empty file). Misidentified reads so far have proven to be low-quality
# or problematic reads, and always get filtered out in subsequent steps, but I
# still like to remove these reads from analyses, but keep them just in case.

### Remove empty sample files --------------------------------------------------
# Make sure all sample files contain reads. Samples with size of 100 bytes or
# below do not have any reads, and this will break the pipeline later if these
# samples are not removed.

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

# This saves the gene1 R1 and R2 fastq sample files only if both the R1 and R2
# sample files have reads.
trimmed_noreads_gene1_R1 <- trimmed_gene1_R1[
  sapply(trimmed_gene1_R1, file.size) < 100
]
file.remove(trimmed_noreads_gene1_R1)
trimmed_noreads_gene1_R2 <- trimmed_gene1_R2[
  sapply(trimmed_gene1_R2, file.size) < 100
]
file.remove(trimmed_noreads_gene1_R2)

print(paste(
  "Here are the samples files for",
  gene1,
  "which contain no reads after primer trimming",
  sep = " "
))
trimmed_noreads_gene1_R1
# This saves the gene2 R1 and R2 fastq sample files only if both the R1 and R2
# sample files have reads.
trimmed_noreads_gene2_R1 <- trimmed_gene2_R1[
  sapply(trimmed_gene2_R1, file.size) < 100
]
file.remove(trimmed_noreads_gene2_R1)
trimmed_noreads_gene2_R2 <- trimmed_gene2_R2[
  sapply(trimmed_gene2_R2, file.size) < 100
]
file.remove(trimmed_noreads_gene2_R2)

print(paste(
  "Here are the samples files for",
  gene2,
  "which contain no reads after primer trimming",
  sep = " "
))
trimmed_noreads_gene2_R1
### Remove mistmatched reads----------------------------------------------------
# Check to see how many wrong-gene occurances there are for gene1. Replace
# "gene1" with your first gene name, and "gene2" with your second gene name
# for all instances below and save the names of the samples
# with these misidentifications.
mismatches_gene1 <- sort(
  list.files(
    trimmed_gene1,
    pattern = gene2,
    full.names = TRUE
  )
)
#mismatches.gene1

# Check to see how many items are in mismatches.gene1
print(paste(
  "Here are the number of reads from which the",
  gene2,
  "primer was removed from samples that were supposed to contain only",
  gene1,
  "amplicons",
  sep = " "
))
length(mismatches_gene1)

# Check the file size of these files to get an estimate of the number of reads
# each micro-contaminate has. If file sizes are < 1kb, it contains less than
# 20 reads (and file sizes below 100 are empty). If you have any files that are
# signficantly larger, you may have contamination issues.
print(paste(
  "Here are the file sizes for trimmed reads from which the",
  gene2,
  "primer was removed from samples that were supposed to contain only",
  gene1,
  "amplicons",
  sep = " "
))
file.size(mismatches_gene1)

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
# Check to see how many items are in mismatches.gene1
print(paste(
  "Here are the number of trimmed reads from which the",
  gene1,
  "primer was removed from samples that were supposed to contain only",
  gene2,
  "amplicons",
  sep = " "
))
length(mismatches_gene2)
# Check the file size of these files to get an estimate of the number of reads
# each micro-contaminate has. If file sizes are < 1kb, it contains less than
# 20 reads (and file sizes below 100 are empty). If you have any files that are
# signficantly larger, you may have contamination issues.
print(paste(
  "Here are the file sizes for trimmed reads from which the",
  gene1,
  "primer was removed from samples that were supposed to contain only",
  gene2,
  "amplicons",
  sep = " "
))
file.size(mismatches_gene2)

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
# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs_gene1 <- sort(
  list.files(
    trimmed_gene1,
    pattern = "_R1.fastq",
    full.names = TRUE
  )
)

fnRs_gene1 <- sort(
  list.files(
    trimmed_gene1,
    pattern = "_R2.fastq",
    full.names = TRUE
  )
)

sample_names_gene1 <- sapply(strsplit(fnFs_gene1, "_trimmed"), `[`, 1)
#fnFs_gene1
#fnRs_gene1
#sample_names_gene1
# Make sure you have the correct number of samples, and that they match the
# number of sample names in the list you made previously.
#length(fnFs_gene1)
#length(fnRs_gene1)
#length(sample_names_gene1)
nsamples_gene1 <- length(sample_names_gene1)
print(paste(
  "Here are the number of",
  gene1,
  "samples that will be analyzed with DADA2:",
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

qualplotF_gene1 <- plotQualityProfile(
  fnFs_gene1[1:nsamples_gene1],
  aggregate = TRUE
)
qualplotF_gene1 <- qualplotF_gene1 +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

# Examine the reverse reads as you did the forward.
qualplotR_gene1 <- plotQualityProfile(
  fnRs_gene1[1:nsamples_gene1],
  aggregate = TRUE
)
qualplotR_gene1 <- qualplotR_gene1 +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

### Export Quality Plots -------------------------------------------------------

ggsave(
  paste0(
    "../data/results/",
    gene1,
    "/qualplotF_",
    gene1,
    ".pdf"
  ),
  plot = qualplotF_gene1,
  width = 9,
  height = 9
)

ggsave(
  paste0(
    "../data/results/",
    gene1,
    "/qualplotR_",
    gene1,
    ".pdf"
  ),
  plot = qualplotR_gene1,
  width = 9,
  height = 9
)

## Gene2 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs_gene2 <- sort(
  list.files(
    trimmed_gene2,
    pattern = "_R1.fastq",
    full.names = TRUE
  )
)
fnRs_gene2 <- sort(
  list.files(
    trimmed_gene2,
    pattern = "_R2.fastq",
    full.names = TRUE
  )
)
sample_names_gene2 <- sapply(strsplit(fnFs_gene2, "_trimmed"), `[`, 1)

# Make sure you have the correct number of samples, and that they match the
# number of sample names in the list you made previously.
#length(fnFs_gene2)
#length(fnRs_gene2)
#length(sample_names_gene2)
nsamples_gene2 <- length(sample_names_gene2)

print(paste(
  "Here are the number of",
  gene1,
  "samples that will be analyzed with DADA2:",
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
qualplotF_gene2 <- plotQualityProfile(
  fnFs_gene2[1:nsamples_gene2],
  aggregate = TRUE
)
qualplotF_gene2 <- qualplotF_gene2 +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

# Examine the reverse reads as you did the forward.
qualplotR_gene2 <- plotQualityProfile(
  fnRs_gene2[1:nsamples_gene2],
  aggregate = TRUE
)
qualplotR_gene2 <- qualplotR_gene2 +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )

### Export Quality Plots -------------------------------------------------------

ggsave(
  paste0(
    "../data/results/",
    gene2,
    "/qualplotF_",
    gene2,
    ".pdf"
  ),
  plot = qualplotF_gene2,
  width = 9,
  height = 9
)

ggsave(
  paste0(
    "../data/results/",
    gene2,
    "/qualplotR_",
    gene2,
    ".pdf"
  ),
  plot = qualplotR_gene2,
  width = 9,
  height = 9
)
