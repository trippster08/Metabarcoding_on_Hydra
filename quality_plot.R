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

## File Housekeeping ===========================================================

# Set up your working directory. If you created your new project in the
# directory you want as your working directory (or came directory from the
# previous step in the pipeline), you don't need to do this, and
# skip to the next RStudio command. If you need to set your working directory,
# substitute your own path for the one below.

# Set a path to the directory with the cutadapt-trimmed reads.
trimmed <- "../data/working/trimmed_reads"
data <- Sys.getenv("data")
numcores <- Sys.getenv("NSLOTS")

# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs <- sort(list.files(trimmed, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(trimmed, pattern = "_R2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Make sure you have the correct number of samples, and that they match the
# number of sample names in the list you made previously.
length(fnFs)
length(fnRs)
length(sample.names)
nsamples <- length(sample.names)

# Make sure all sample files contain reads. Samples with size of 50 bytes or
# below do not have any reads, and this will break the pipeline later if these
# samples are not removed.
file.size(fnFs)

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnFs.exists <- fnFs[file.size(fnFs) > 50 & file.size(fnRs) > 50]
length(fnFs.exists)

# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnRs.exists <- fnRs[file.size(fnFs) > 50 & file.size(fnRs) > 50]
length(fnRs.exists)
file.size(fnFs.exists)

# Redefine fnFs and fnRs as only the existing read files, and check
fnFs <- fnFs.exists
fnRs <- fnRs.exists
length(fnFs)
length(fnRs)
file.size(fnFs)

# Update your samples names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
nsamples <- length(sample.names)
head(sample.names)

## Filter and Trim =============================================================

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
qualplotF <- plotQualityProfile(
  fnFs[1:nsamples],
  aggregate = TRUE
)
qualplotF <- qualplotF + scale_x_continuous(limits = c(100, 300), breaks = seq(100, 300, 10))
ggsave(
  "../data/working/qualplotF.pdf",
  plot = qualplotF,
  width = 9,
  height = 9
)


# Examine the reverse reads as you did the forward.
qualplotR <- plotQualityProfile(
  fnRs[1:nsamples],
  aggregate = TRUE
)
qualplotR <- qualplotR + scale_x_continuous(limits = c(100, 300), breaks = seq(100, 300, 10))
ggsave(
  "../data/working/qualplotR.pdf",
  plot = qualplotR,
  width = 9,
  height = 9
)