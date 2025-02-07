# TRIM FILTER DENOISE MERGE ####################################################
library(dada2)
library(tidyverse)
library(seqinr)
library(digest)

## Trim Reads ==================================================================
numcores <- Sys.getenv("NSLOTS")
args <- commandArgs(trailingOnly = TRUE)
trimmed <- args[1]
truncF <- as.numeric(args[2])
truncR <- as.numeric(args[3])


# This creates two vectors. One contains the names for forward reads (R1, called
# fnFs) and the other for reverse reads (R2, called fnRs).
fnFs <- sort(list.files(trimmed, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(trimmed, pattern = "_R2.fastq.gz", full.names = TRUE))

# Create a list of sample names
sample.names <- sapply(strsplit(basename(fnFs), "_trimmed"), `[`, 1)
print("Here are the sample names of the first 6 trimmed R1 files:")

head(sample.names)

# If you have sample files with no reads, you must remove both the forward and
# reverse reads, regardless if one has reads (although if one is empty,
# the other should be as well).

### Remove empty sample files --------------------------------------------------
# This saves the R1 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnFs.exists <- fnFs[file.size(fnFs) > 100 & file.size(fnRs) > 100]


# This saves the R2 fastq for the sample file only if both the R1 and R2 sample
# files have reads.
fnRs.exists <- fnRs[file.size(fnFs) > 100 & file.size(fnRs) > 100]


# Redefine fnFs and fnRs as only the existing read files, and check
fnFs <- fnFs.exists
fnRs <- fnRs.exists
file.size(fnFs)

# Update your samples names
sample.names <- sapply(strsplit(basename(fnFs), "_trimmed"), `[`, 1)

# This creates files for the reads that will be quality filtered with dada2
# in the next step.
filtFs <- file.path(
  "../data/working",
  "filtered",
  paste0(
    sample.names,
    "_F_filt.fastq.gz"
  )
)
filtRs <- file.path(
  "../data/working",
  "filtered",
  paste0(
    sample.names,
    "_R_filt.fastq.gz"
  )
)

# This inserts sample names to these newly created files.
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
  multithread = TRUE,
  verbose = TRUE
)

# Usually we don't have that many samples, so I just look at "out" in its
# entirety, but if there are lots of samples, just look at the first 6.
out

# Export out as a tsv
write.table(
  out,
  file = "../data/results/filtered_reads.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

# Save all the objects created to this point
save(
  trimmed,
  truncF,
  truncR,
  fnFs,
  fnRs,
  sample.names,
  filtFs,
  filtRs,
  out,
  file = "../data/results/out.Rdata"
)


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
print("Here are the number of samples that were removed because they no longer contain reads after filtering")
length(filtFs) - length(filtFs[exists])

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
  "../data/results/errorplotsF.pdf",
  plot = error.plots.F,
  width = 9,
  height = 9
)
ggsave(
  "../data/results/errorplotsR.pdf",
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

# Save all the objects created between out and here
save(
  exists,
  filtFs,
  filtRs,
  errF,
  errR,
  dadaFs,
  dadaRs,
  file = "../data/results/denoise.feattab.RData"
)



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

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2. I've added a column to the typical
# tutorial version of this that gives us the percentage of reads that made it
# through the process. This is a good quick way to see if something is wrong
# (i.e. only a small proportion make it through).
getN <- function(x) sum(getUniques(x))
track <- cbind(
  out, 
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(merged, getN),
  rowSums(seqtab.nochim),
  100 * (rowSums(seqtab.nochim) / out[, 1]))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "%kept")
rownames(track) <- sample.names

# Export this table as a .tsv
write.table(
  track,
  file = "../data/results/track_reads.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Examine Sequence Lengths ====================================================

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
seq.length.table <- table(nchar(getSequences(seqtab.nochim)))
# Export this table as a .tsv
write.table(
  seq.length.table,
  file="../data/results/ASV_lengths_table.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Create And Use md5 Hash =====================================================
# To create a Sequence list with md5 hash instead of ASVs, we first need to
# create a list of md5 hash's of all ASV's.

# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2. We are going to use this list to create md5 hashes. Use whatever
#  table you will later use for your analyses (e.g. seqtab.nochim)
repseq <- getSequences(seqtab.nochim)

# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
repseq.md5 <- c()
for (i in seq_along(repseq)) {
  repseq.md5[i] <- digest(
    repseq[i],
    serialize = FALSE,
    algo = "md5"
  )
}

# Add md5 hash to the sequence-table from the DADA2 analysis.
seqtab.nochim.md5 <- seqtab.nochim
colnames(seqtab.nochim.md5) <- repseq.md5

# Create an md5/ASV table, with each row as an ASV and it's representative md5
# hash.
repseq.md5.asv <- tibble(repseq.md5, repseq)
# Rename column headings
colnames(repseq.md5.asv) <- c("md5", "ASV")

# Transpose the sequence-table, and convert the result into a tibble.
seqtab.nochim.transpose.md5 <- as_tibble(t(seqtab.nochim.md5), rownames = "ASV")

# Save all the objects created between denoise and here
save(
  merged,
  seqtab,
  seqtab.nochim,
  getN,
  track,
  seq.length.table,
  repseq,
  repseq.md5,
  seqtab.nochim.md5,
  repseq.md5.asv,
  seqtab.nochim.transpose.md5,
  file = "../data/results/feattab.RData"
)


## Export Feature-Table with md5 Hash =========================================
# This exports a feature-table: row of ASV's (shown as a md5 hash instead
# of sequence), columns of samples, and values = number of reads. With this
# table you will also need a file that relates each ASV to it's representative
# md5 hash. We download this in the next section.

write.table(
  seqtab.nochim.transpose.md5,
  file = "../data/results/feature-table_md5.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Export Representative Sequences table/fasta =================================
# Here we export our our representative sequences, either as a fasta (with the
# md5 hash as the ASV name), or as a table with ASV and md5 hash as columns.

# This exports all the ASVs in fasta format, with ASV hash as the sequence
# name. This is analogous to the representative sequence output in Qiime2.
write.fasta(
  sequences = as.list(repseq.md5.asv$ASV),
  names = repseq.md5.asv$md5,
  open = "w",
  as.string = FALSE,
  file.out = "../data/results/rep-seq.fas"
)

# This exports all the ASVs and their respective md5 hashes as a two-column
# table.
write.table(
  repseq.md5.asv,
  file = "../data/results/representative_sequence_table_md5.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)