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
trimmed <- "../data/working/trimmed_reads"
data <- "../data/"
trimmed_R1 <- sort(
  list.files(
    trimmed,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
trimmed_R2 <- sort(
  list.files(
    trimmed,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
trimmed_noreads_R1 <- trimmed_R1[sapply(trimmed_R1, file_size) < 100]
file.remove(trimmed_noreads_R1)
trimmed_noreads_R2 <- trimmed_R2[sapply(trimmed_R2, file_size) < 100]
file.remove(trimmed_noreads_R2)

print(
  "Here are the samples files for which contain no reads after primer trimming"
)
trimmed_noreads_R1

# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs <- sort(
  list.files(
    trimmed,
    pattern = "_R1.fastq.gz",
    full.names = TRUE
  )
)
fnRs <- sort(
  list.files(
    trimmed,
    pattern = "_R2.fastq.gz",
    full.names = TRUE
  )
)

# Make your samples names
sample_names <- sapply(strsplit(fnFs, "_trimmed"), `[`, 1)
nsamples <- length(sample_names)

## Make Quality Plots ----------------------------------------------------------

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
qualplotF <- qualplotF +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )
ggsave(
  "../data/results/qualplotF.pdf",
  plot = qualplotF,
  width = 9,
  height = 9
)

# Examine the reverse reads as you did the forward.
qualplotR <- plotQualityProfile(
  fnRs[1:nsamples],
  aggregate = TRUE
)
qualplotR <- qualplotR +
  scale_x_continuous(
    limits = c(100, 300),
    breaks = seq(100, 300, 10)
  )
ggsave(
  "../data/results/qualplotR.pdf",
  plot = qualplotR,
  width = 9,
  height = 9
)
