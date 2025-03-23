# TRIM FILTER DENOISE MERGE ####################################################
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## Trimmed Sequences ===========================================================

args <- commandArgs(trailingOnly = TRUE)
path_to_trimmed <- args[1]
truncF <- as.numeric(args[2])
truncR <- as.numeric(args[3])

# Save project name as an object
project_name <- basename(dirname(getwd()))
print(paste0("This project is named ", project_name))

# Load the RData from "quality_plot_multigene.R"
load("../data/working/1_trim_qual.RData")

# This creates a vector of the path for forward reads (R1, called trimmed_F).
trimmed_F <- sort(list.files(
  path_to_trimmed,
  pattern = "_R1.fastq.gz",
  full.names = TRUE
))

# Make a vector of sample names from your trimmed reads.
sample_names_trimmed <- sapply(
  strsplit(basename(trimmed_F), "_S\\d{1,3}_"),
  `[`,
  1
)

# Give the vectors names
names(trimmed_F) <- sample_names_trimmed
names(trimmed_R) <- sample_names_trimmed


## Get Read Counts of Trimmed Samples ==========================================
# Count the number of reads in each trimmed sample. Since cutadapt only
# keeps paired reads, we only need to count forward samples.
sequence_counts_trimmed <- sapply(trimmed_F, function(file) {
  fastq_data <- readFastq(file)
  length(fastq_data)
})
names(sequence_counts_trimmed) <- sample_names_trimmed

# This creates files for the reads that will be quality filtered with dada2
# in the next step.
filtered_F <- file.path(
  "../data/working",
  "filtered_sequences",
  paste0(
    sample_names_trimmed,
    "_F_filt.fastq.gz"
  )
)
filtered_R <- file.path(
  "../data/working",
  "filtered_sequences",
  paste0(
    sample_names_trimmed,
    "_R_filt.fastq.gz"
  )
)

# This inserts sample names to these newly created files.
names(filtered_F) <- sample_names_trimmed
names(filtered_R) <- sample_names_trimmed

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(X,Y)" is how you tell Dada2 where to truncate all forward (X) and
# reverse (Y) reads. Using "0" means reads will not be truncated.
# maxEE sets how many expected errors are allowed before a read is filtered out.

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 4,4) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.

filtered_summary <- filterAndTrim(
  trimmed_F,
  filtered_F,
  trimmed_R,
  filtered_R,
  truncLen = c(truncF, truncR),
  maxN = 0,
  maxEE = c(4, 4),
  rm.phix = TRUE,
  truncQ = 2,
  compress = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

# After filtering, if you have any samples with no reads, you much remove them.
# This step changes filtered_F and filtered_R to only contain the names of
# samples with reads. Do this only if there are samples in
# "filtered_summary" with zero reads.
exists <- file.exists(filtered_F) & file.exists(filtered_R)
print(
  "Here are the number of samples that were removed because they no longer contain reads after filtering"
)
length(filtered_F) - length(filtered_F[exists])

filtered_F <- filtered_F[exists]
filtered_R <- filtered_R[exists]


# Set a path to the directory with the dada2-filtered reads.
path_to_filtered <- "../data/working/filtered_sequences"

# Get sample names for filtered reads
sample_names_filtered <- sapply(
  strsplit(basename(filtered_F), "_[FR]_filt"),
  `[`,
  1
)

# Count how many reads remain in each sample after filtering
sequence_counts_filtered <- sapply(filtered_F, function(file) {
  fastq_data <- readFastq(file)
  length(fastq_data)
})

# Name the counts with sample names
names(sequence_counts_filtered) <- sample_names_filtered

print("Here are the read counts for each filtered sample:")
sequence_counts_filtered

# Save all the objects created to this point
save(
  path_to_trimmed,
  truncF,
  truncR,
  trimmed_F,
  sample_names_trimmed,
  filtered_F,
  filtered_R,
  filtered_summary,
  path_to_filtered,
  sample_names_filtered,
  sequence_counts_filtered,
  file = "../data/working/filtered_summary.RData"
)

# Export filtered_summary as a tsv
write.table(
  filtered_summary,
  file = paste0("../data/results/", project_name, "_filtered_read_count.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Error Estimation ============================================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errors_F <- learnErrors(
  filtered_F,
  nbases = 1e+08,
  errorEstimationFunction = loessErrfun,
  multithread = TRUE,
  randomize = FALSE,
  MAX_CONSIST = 10,
  OMEGA_C = 0,
  qualityType = "Auto",
  verbose = FALSE
)

errors_R <- learnErrors(
  filtered_R,
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
error_plots_F <- plotErrors(errors_F, nominalQ = TRUE)
error_plots_R <- plotErrors(errors_R, nominalQ = TRUE)

ggsave(
  paste0("../data/results/", project_name, "_error_plots_F.pdf"),
  plot = error_plots_F,
  width = 9,
  height = 9
)
ggsave(
  paste0("../data/results/", project_name, "_error_plots_R.pdf"),
  plot = error_plots_R,
  width = 9,
  height = 9
)

# Save the objects created since filtered_summary
save(
  errors_F,
  errors_R,
  error_plots_F,
  error_plots_R,
  file = "../data/working/errors.RData"
)
## Denoising ===================================================================

# This applies the "core sample inference algorithm" (i.e. denoising) in dada2
# to get corrected unique sequences. The two main inputs are the first, which is
# the filtered sequences (filtered_F), and "err =" which is the error file from
# learnErrors (effF).
denoised_F <- dada(
  filtered_F,
  err = errors_F,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

denoised_R <- dada(
  filtered_R,
  err = errors_R,
  errorEstimationFunction = loessErrfun,
  selfConsist = FALSE,
  pool = FALSE,
  multithread = TRUE,
  verbose = TRUE
)

# Save the denoise objects
save(
  denoised_F,
  denoised_R,
  file = "../data/working/denoise.RData"
)

## Merge Paired Sequences ======================================================

# Here we merge the paired reads. merged calls for the forward denoising result
# (denoised_F), then the forward filtered and truncated reads (filtered_F),
# then the same for the reverse reads (denoised_R and filtered_R).

# You can change the minimum overlap (minOverlap), and the number of mismatches
# that are allowed in the overlap region (maxMismatch). Default values are
# shown.

# mergePairs results in a data.frame from each sample that contains a row for
# "...each unique pairing of forward/reverse denoised sequences." The data.frame
# also contains multiple columns describing data for each unique merged
# sequence.
merged_reads <- mergePairs(
  denoised_F,
  filtered_F,
  denoised_R,
  filtered_R,
  minOverlap = 12,
  maxMismatch = 0,
  verbose = TRUE
)

## Create Sequence-Table =======================================================

# Now we make a sequence-table containing columns for unique sequence (ASV),
# rows for each sample, and numbers of reads for each ASV/sample combination as
# the body of the table. This table is in the reverse orientation of the feature
# table output by Qiime2 (which uses samples as columns), but we easily
# transpose this table if needed later (and we will later). I think for now, I
# will use "sequence-table" for the table with columns of sequences, and
# "feature-table" for tables with columns of samples.
seqtab <- makeSequenceTable(merged_reads)
# This describes the dimensions of the table just made
print(paste(
  "These are the dimensions of your newly created Sequence-Table:",
  dim(seqtab),
  sep = " "
))


## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

print(paste(
  "These are the dimensions of your chimera-free Sequence-Table:",
  dim(seqtab),
  sep = " "
))

# Make a list of the ASVs that are considered chimeras, in case you want to look
# at them later
chimeras_list <- isBimeraDenovoTable(
  seqtab,
  multithread = TRUE,
  verbose = TRUE
)
# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2.
repseq_all <- getSequences(seqtab)
# Get a list of just chimera ASVs
repseq_chimera <- repseq_all[chimeras_list]
# Make and add md5 hash to the repseq_chimera
repseq_chimera_md5 <- c()
for (i in seq_along(repseq_chimera)) {
  repseq_chimera_md5[i] <- digest(
    repseq_chimera[i],
    serialize = FALSE,
    algo = "md5"
  )
}

# Export this as a fasta
write.fasta(
  sequences = as.list(repseq.chimera),
  names = repseq.chimera_md5,
  open = "w",
  as.string = FALSE,
  file.out = paste0("../data/results/", project_name, "_rep-seq_chimeras.fas")
)

## Examine Sequence Lengths ====================================================

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.
seq_length_table <- table(nchar(getSequences(seqtab_nochim)))
# Export this table as a .tsv
write.table(
  seq_length_table,
  file = paste0("../data/results/", project_name, "_ASV_lengths_table.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2 or "track" from the DADA2 tutorial.
# This is a good quick way to see if something is wrong (i.e. only a small
# proportion make it through).

# First make a table for the post-filtered samples, including denoised,
# merged, and non-chimeric read counts
getN <- function(x) sum(getUniques(x))
sequence_counts_postfiltered <- as_tibble(
  cbind(
    sapply(denoised_F, getN),
    sapply(denoised_R, getN),
    sapply(merged_reads, getN),
    rowSums(seqtab_nochim)
  ),
  .name_repair = "unique"
) %>%
  mutate(Sample_ID = sample_names_filtered) %>%
  select(
    Sample_ID,
    Denoised_Reads_F = ...1,
    Denoised_Reads_R = ...2,
    Merged_Reads = ...3,
    Non_Chimeras = ...4
  )

# Then we are going to add the postfiltered read count data to the three count
# data objects we already have (raw, trimmed, filtered).
track_reads <- tibble(
  sequence_counts_raw,
  Sample_ID = names(sequence_counts_raw)
) %>%
  left_join(
    tibble(sequence_counts_trimmed, Sample_ID = names(sequence_counts_trimmed)),
    join_by(Sample_ID)
  ) %>%
  left_join(
    tibble(
      sequence_counts_filtered,
      Sample_ID = names(sequence_counts_filtered)
    ),
    join_by(Sample_ID)
  ) %>%
  left_join(
    sequence_counts_postfiltered,
    join_by(Sample_ID)
  ) %>%
  mutate(Proportion_Kept = Non_Chimeras / sequence_counts_raw) %>%
  select(Sample_ID, everything()) %>%
  select(
    Sample_ID,
    Raw_Reads = sequence_counts_raw,
    Trimmed_Reads = sequence_counts_trimmed,
    Filtered_Reads = sequence_counts_filtered,
    everything()
  )


# Export this table as a .tsv
write.table(
  track_reads,
  file = paste0("../data/results/", project_name, "_track_reads.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

## Export Sequence-Table =======================================================
# This exports a sequence-table: columns of ASV's, rows of samples, and
# values = number of reads. This is the only export you need for downstream
# analyses. You can do anything you want in this pipeline with this table. This
# table can also be easily merged with other tables from the same project but
# from different runs before downstream analyses.

# If you have mulitple Miseqruns for the same project that will need to be
# combined for further analyses, you may want to name this file
# "PROJECTNAME_MISEQRUN_sequence-table.tsv" to differentiate different runs.
# In "5 Metabarcoding_R_Pipeline_RStudio_ImportCombine" we'll show how to
# combine data from separate runs for analyses.

write.table(
  seqtab_nochim,
  file = paste0("../data/results/", project_name, "_sequence-table.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)
# NOTE!!!
# The Sequence-Table in this format is very unwieldy, since each column name is
# an entire ASV. Instead, we can convert each ASV into a short "hash" using
# the md5 encryption model, creating a 32bit representative of each ASV. Each
# hash is essentially unique to the ASV it is representing. We would then
# replace the ASVs in the column headings with their representative md5 hash.
# However, having an ASV hash as a column heading requires the creation of a
# Representative Sequence list, which tells us which hash goes with which ASV.
# gives the user a representative-sequence fasta that contains the ASV, labelled
# with its specfic md5 hash.
# If you want to export a Sequence-Table with a md5 hash instead of ASV sequence
# for each ASV, skip this and go to the next section.

## Create And Use md5 Hash =====================================================
# To create a Sequence list with md5 hash instead of ASVs, we first need to
# create a list of md5 hash's of all ASV's.

# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2. We are going to use this list to create md5 hashes. Use whatever
#  table you will later use for your analyses (e.g. seqtab_nochim)
repseq_nochim <- getSequences(seqtab_nochim)

# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
repseq_nochim_md5 <- c()
for (i in seq_along(repseq_nochim)) {
  repseq_nochim_md5[i] <- digest(
    repseq_nochim[i],
    serialize = FALSE,
    algo = "md5"
  )
}

# Add md5 hash to the sequence-table from the DADA2 analysis.
seqtab_nochim_md5 <- seqtab_nochim
colnames(seqtab_nochim_md5) <- repseq_nochim_md5

# Export this sequence table with column headings as md5 hashs instead of ASV
# sequences
write.table(
  seqtab_nochim_md5,
  file = paste0("../data/results/", project_name, "_sequence-table-md5.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

# Create an md5/ASV table, with each row as an ASV and it's representative md5
# hash.
repseq_nochim_md5_asv <- tibble(repseq_nochim_md5, repseq)
# Rename column headings
colnames(repseq_nochim_md5_asv) <- c("md5", "ASV")


## Create and Export Feature-Table =============================================
# This creates and exports a feature-table: row of ASV's (shown as a md5 hash
# instead of sequence), columns of samples, and values = number of reads. With
# this table you will also need a file that relates each ASV to it's
# representative md5 hash. We download this in the next section.

# Transpose the sequence-table, and convert the result into a tibble.
seqtab_nochim_transpose_md5 <- as_tibble(
  t(seqtab_nochim_md5),
  rownames = "ASV"
)


write.table(
  seqtab_nochim_transpose_md5,
  file = paste0("../data/results/", project_name, "_feature-table_md5.tsv"),
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
  sequences = as.list(repseq_nochim_md5_asv$ASV),
  names = repseq_nochim_md5_asv$md5,
  open = "w",
  as.string = FALSE,
  file.out = paste0("data/results/", project_name, "_rep-seq.fas")
)

# This exports all the ASVs and their respective md5 hashes as a two-column
# table.
write.table(
  repseq_nochim_md5_asv,
  file = paste0(
    "data/results/",
    project_name,
    "_representative_sequence_md5_table.tsv"
  ),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# Save all the objects created between denoise and here
save(
  merged_reads,
  seqtab,
  seqtab_nochim,
  track_reads,
  chimeras_list,
  repseq_all,
  repseq_chimera,
  seq_length_table,
  repseq_nochim,
  repseq_nochim_md5,
  seqtab_nochim_md5,
  repseq_nochim_md5_asv,
  seqtab_nochim_transpose_md5,
  file = "../data/working/feattab.RData"
)

## Create a Sequence-List Table ================================================
# This creates a table containing three columns: sample name, ASV, and read
# count. Each row is a separate sample/ASV combination. This is a tidy table,
# which means each column contains a single variable and each rowh a single
# observation. This is a good table format for storage of DADA2 results because
# it can be easily concatenated with other sequence-list tables in Excel or any
# text-editing software (unlike the sequence-table), yet it still contains all
# the information needed from our trimming and denoising steps.

# You will also need to make a Sequence-List table if you want to export a
# feature-to-fasta file later. A feature-to-fasta file contains every
# combination of ASV and sample, and each sequence is named with the sample
# name, md5 hash and number of reads of that ASV in that sample. It's a good
# way to look at ASV distributions in a phylogenetic tree.

# Convert the sequence-table from your DADA2 results into a tibble,
# converting row names to column 1, labeled "sample". A tibble is a more
# versatile data.frame, but it does not have row headings
# (among other differences, see https://tibble.tidyverse.org/). We'll need
# this to be a tibble for the next step.
seqtab_nochim_tb <- as_tibble(seqtab_nochim, rownames = "sample")

# The sequence-table has a column with sample names, and N columns of ASV's
# containing count values. We want all the count data to be in a single column,
# so we use a tidyr command called "pivot_longer" to make the table "tall",
# which means the table goes from 17 x 2811 to 47770 x 3 for example
# (47770 = 2810 x 17. 2810 instead of 2811 because the first column of the
# original table contains sample names, not counts). This makes the table tidier
# (meaning that each column is now a true variable).
seqtab_nochim_tall <- seqtab_nochim_tb %>%
  pivot_longer(
    !sample,
    names_to = "ASV",
    values_to = "count"
  )

# Remove rows with sequence counts = 0. This removes any sample in which a
# particular ASV was not found.
seqtab_nochim_tall_nozero <- subset(seqtab_nochim_tall, count != 0)

## Create and Export feature-to-fasta ==========================================
# This creates a fasta file containing all the ASV's for each sample. Each ASV
# will be labeled with the sample name, ASV hash, and number of reads of that
# ASV in that sample. This was derived from a python script from Matt Kweskin
# called featuretofasta.py (hence the name).

# Save the ASV sequences from the sequence-list table
# (seqtab_nochim_tall_nozera) as a new list.
repseq_tall <- seqtab_nochim_tall_nozero$ASV

# Convert the sequences into md5 hashs, as we did earlier. md5 hashs are
# consistent across jobs, meaning identical sequences from different projects or
# being converted by different programs will result in the same hash (i.e.
# hashs here will match hashs above)
repseq_tall_md5 <- c()
for (i in seq_along(repseq_tall)) {
  repseq_tall_md5[i] <- digest(
    repseq_tall[i],
    serialize = FALSE,
    algo = "md5"
  )
}

# Attach the ASV hashes as a column (called "md5") to the tall table. The
# table should now have 4 columns, and each row of the "md5" column should
# be a md5 hash of its respective ASV.
seqtab_nochim_tall_nozero_md5 <- cbind(
  seqtab_nochim_tall_nozero,
  md5 = repseq_tall_md5
)

# Create a new column in this table that contains "sample", "feature", and
# "count", concatenated. This is the heading for each sequence in the fasta file
# created by Matt Kweskin's script "featuretofasta.py"
seqtab_nochim_tall_nozero_md5_ftf <- seqtab_nochim_tall_nozero_md5 %>%
  mutate(sample_md5_count = paste(sample, md5, count, sep = "_")) %>%
  select(sample, md5, count, sample_md5_count, ASV)

# Create a fasta-formatted file of each row sequence (i.e. ASV), with a heading
# of "sample_feature_count".
write.fasta(
  sequences = as.list(seqtab_nochim_tall_nozero_md5_ftf$ASV),
  names = seqtab_nochim_tall_nozero_md5_ftf$sample_feature_count,
  open = "w",
  as.string = FALSE,
  file.out = paste0("data/results/", project_name, "_feature-to-fasta.fas")
)

save.image(file = "../data/working/2_filter_denoise_merge.RData")
