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
numcores <- Sys.getenv("NSLOTS")
gene1 <- args[1]
gene2 <- args[2]
trimmed.gene1 <- paste0("../data/working/trimmed_reads/", gene1)
trimmed.gene2 <- paste0("../data/working/trimmed_reads/", gene2)

# We will now run the rest of this section in multiple times, once for each
# gene present in the Illumina run. Again, replace each "gene1", "gene2",
# "gene3", etc with your specific gene name.

## Gene1 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs.gene1 <- sort(
  list.files(
    trimmed.gene1,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
fnRs.gene1 <- sort(
  list.files(
    trimmed.gene1,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)
sample.names.gene1 <- sapply(strsplit(fnFs.gene1, "_"), `[`, 1)

# Make sure you have the correct number of samples, and that they match the
# number of sample names in the list you made previously.
length(fnFs.gene1)
length(fnRs.gene1)
length(sample.names.gene1)
nsamples.gene1 <- length(sample.names.gene1)

# Make sure all sample files contain reads. Samples with size of 50 bytes or
# below do not have any reads, and this will break the pipeline later if these
# samples are not removed.
file.size(fnFs.gene1)

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnFs.exists.gene1 <- fnFs.gene1[
  file.size(fnFs.gene1) > 50 & file.size(fnRs.gene1) > 50
]
length(fnFs.exists.gene1)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnRs.exists.gene1 <- fnRs.gene1[
  file.size(fnFs.gene1) > 50 & file.size(fnRs.gene1) > 50
]
length(fnRs.exists.gene1)
file.size(fnFs.exists.gene1)

# Redefine fnFs and fnRs as only the existing read files, and check
fnFs.gene1 <- fnFs.exists.gene1
fnRs.gene1 <- fnRs.exists.gene1
length(fnFs.gene1)
length(fnRs.gene1)
file.size(fnFs.gene1)

# Update your samples names
sample.names.gene1 <- sapply(strsplit(basename(fnFs.gene1), "_"), `[`, 1)
nsamples.gene1 <- length(sample.names.gene1)
head(sample.names.gene1)

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
qualplotF.gene1 <- plotQualityProfile(
  fnFs.gene1[1:nsamples.gene1],
  aggregate = TRUE
)
qualplotF.gene1 <- qualplotF.gene1 +
scale_x_continuous(
  limits = c(100, 300),
  breaks = seq(100, 300, 10)
)

# Examine the reverse reads as you did the forward.
qualplotR.gene1 <- plotQualityProfile(
  fnRs.gene1[1:nsamples.gene1],
  aggregate = TRUE
)
qualplotR.gene1 <- qualplotR.gene1 +
scale_x_continuous(
  limits = c(100, 300), 
  breaks = seq(100, 300, 10)
)

### Export Quality Plots -------------------------------------------------------

ggsave(
  paste0("../data/working/qualplotF_", gene1, ".pdf"),
  plot = qualplotF.gene1,
  width = 9,
  height = 9
)

ggsave(
  paste0("../data/working/qualplotR_", gene1, ".pdf"),
  plot = qualplotR.gene1,
  width = 9,
  height = 9
)

## Gene2 =======================================================================
# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs.gene2 <- sort(
  list.files(
    trimmed.gene2,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
fnRs.gene2 <- sort(
  list.files(
    trimmed.gene2,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)
sample.names.gene2 <- sapply(strsplit(fnFs.gene2, "_"), `[`, 1)

# Make sure you have the correct number of samples, and that they match the
# number of sample names in the list you made previously.
length(fnFs.gene2)
length(fnRs.gene2)
length(sample.names.gene2)
nsamples.gene2 <- length(sample.names.gene2)

# Make sure all sample files contain reads. Samples with size of 50 bytes or
# below do not have any reads, and this will break the pipeline later if these
# samples are not removed.
file.size(fnFs.gene2)

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnFs.exists.gene2 <- fnFs.gene2[
  file.size(fnFs.gene2) > 50 & file.size(fnRs.gene2) > 50
]
length(fnFs.exists.gene2)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnRs.exists.gene2 <- fnRs.gene2[
  file.size(fnFs.gene2) > 50 & file.size(fnRs.gene2) > 50
]
length(fnRs.exists.gene2)
file.size(fnFs.exists.gene2)

# Redefine fnFs and fnRs as only the existing read files, and check
fnFs.gene2 <- fnFs.exists.gene2
fnRs.gene2 <- fnRs.exists.gene2
length(fnFs.gene2)
length(fnRs.gene2)
file.size(fnFs.gene2)

# Update your samples names
sample.names.gene2 <- sapply(strsplit(basename(fnFs.gene2), "_"), `[`, 1)
nsamples.gene2 <- length(sample.names.gene2)
head(sample.names.gene2)

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
qualplotF.gene2 <- plotQualityProfile(
  fnFs.gene2[1:nsamples.gene2],
  aggregate = TRUE
)
qualplotF.gene2 <- qualplotF.gene2 +
scale_x_continuous(
  limits = c(100, 300),
  breaks = seq(100, 300, 10)
)

# Examine the reverse reads as you did the forward.
qualplotR.gene2 <- plotQualityProfile(
  fnRs.gene2[1:nsamples.gene2],
  aggregate = TRUE
)
qualplotR.gene2 <- qualplotR.gene2 +
scale_x_continuous(
  limits = c(100, 300),
  breaks = seq(100, 300, 10)
)

### Export Quality Plots -------------------------------------------------------

ggsave(
  paste0("../data/working/qualplotF_", gene2, ".pdf"),
  plot = qualplotF.gene2,
  width = 9,
  height = 9
)

ggsave(
  paste0("../data/working/qualplotR_", gene2, ".pdf"),
  plot = qualplotR.gene2,
  width = 9,
  height = 9
)