# TRIM FILTER DENOISE MERGE ####################################################
library(dada2)
library(tidyverse)
library(seqinr)

## Trim Reads ==================================================================

trimmed <- "data/working/trimmed_reads"
data <- Sys.getenv("data")
truncF <- Sys.getenv("truncF")
truncR <- Sys.getenv("truncR")
numcores <- Sys.getenv("NSLOTS")

# Create a list of sample names
trimmed.reads <- list.files(trimmed[str_detect(trimmed, "R1.fastq.gz")])
sample.names <- sapply(strsplit(basename(trimmed.reads), "_"), `[`, 1)


# This creates files for the reads that will be quality filtered with dada2
# in the next step.
filtFs <- file.path(trimmed, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(trimmed, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# This inserts sample names to these newly created files. You'll notice that in
# the environment pane, the description of filtFs and filtRs goes from
# "chr [1:N]" to "Named chr [1:N]"
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(i,j)" is how you tell Dada2 where to truncate all forward (i) and
# reverse (j) reads. Using "0" means reads will not be truncated.
# maxEE sets how many expected errors are allowed before a read is filtered out.

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 6,6) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.

out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  truncLen = c(truncF, truncR),
  maxN = 0,
  maxEE = c(6, 6),
  rm.phix = TRUE,
  truncQ = 2,
  compress = TRUE,
  multithread = 8
)

# Usually we don't have that many samples, so I just look at "out" in its
# entirety, but if there are lots of samples, just look at the first 6.
out
head(out)

# After filtering, if there are any samples that have no remaining reads
# (i.e. reads.out = 0), you will get the following error running learnErrors:
# "Error in derepFastq(fls[[i]], qualityType = qualityType) : Not all provided
# files exist." That is because while these empty sample names still exist in
# filtFs and filtRs, there is no data connected to these names, and dada2
# doesn't like that.

# This step changes filtFs and filtRs to only contain the names of samples with
# reads. Do this only if there are samples in "out" with zero reads.

# You will notice that the number of items in filtFs is now the number of
# samples with reads (i.e. the description for filtFs and filtRs goes from
# "Named chr [1:N]" to "Named chr [1:N-(# of empty samples)]).
exists <- file.exists(filtFs) & file.exists(filtRs)
length(exists)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]


## Error Estimation ============================================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errF <- learnErrors(
  filtFs,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

errR <- learnErrors(
  filtRs,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

# We can visualize the estimated error rates to make sure they don't look too
# crazy. The red lines are error rates expected under the "...nominal defintion
# of the Q-score." The black dots are "...observed error rates for each
# consensus quality score." The black line shows the "...estimated error rates
# after convergence of the machine-learning algorithm." I think the main things
# to look at here are to make sure that each black line is a good fit to the
# observed error rates, and that estimated error rates decrease with increased
# quality.
error.plots.F <- plotErrors(errF, nominalQ = TRUE)
error.plots.R <- plotErrors(errR, nominalQ = TRUE)

ggsave(
  "../data/working/errorplotsF.pdf",
  plot = error.plots.F,
  width = 9,
  height = 9
)
ggsave(
  "../data/working/errorplotsR.pdf",
  plot = error.plots.R,
  width = 9,
  height = 9
)

## Denoising ===================================================================

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtFs), and "err =" which is the error file from
# learnErrors (effF).
dadaFs <- dada(
  filtFs,
  err = errF,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

dadaRs <- dada(
  filtRs,
  err = errR,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# This looks at the dada-class list of objects that was created by the "dada"
# command. It gives a brief summary of the denoising results, and gives some
# parameters values used.
dadaFs[[1]]
dadaRs[[1]]

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. merged calls for the forward denoising result
# (dadaFs), then the forward filtered and truncated reads (filtFs), then the
# same for the reverse reads (dadaRs and filtRs).

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.
merged <- mergePairs(
  dadaFs,
  filtFs,
  dadaRs,
  filtRs,
  minOverlap = 12,
  maxMismatch = 0,
  verbose = TRUE
)

# Inspect the merged sequences from the data.frame of the first sample (and the
# 6th sample).
head(merged[[1]])
head(merged[[6]])

## Create Feature-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab <- makeSequenceTable(merged)
# This describes the dimensions of the table just made
dim(seqtab)

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
table(nchar(getSequences(seqtab)))

# If we want to remove the "excessively" long or short sequences, we can do so
# here. The lengths used here are arbitrary. I'm not sure how to justify a
# cut-off, to be honest. You can sometimes see the a pattern here corresponding
# to codon position for protein-coding genes (there are more ASV's in multiples
# of three), so you may cut where the pattern is no longer visible (i.e. there
# are not more reads in lengths in multiples of threes than at other lengths).
# I tend not to remove any ASV's at this point

# In this example, we only keep reads between 298 and 322 bp in length.
seqtab313 <- seqtab[, nchar(colnames(seqtab)) %in% 298:322]
dim(seqtab313)
table(nchar(getSequences(seqtab313)))

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)
# We look at the dimensions of the new sequence-table
dim(seqtab.nochim)