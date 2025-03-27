# FILTER DENOISE MERGE #########################################################
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## Track Reads Through Dada2 Process ===========================================

# Load the RData from "quality_plot_multigene.R"
load("../data/working/7_chimera.RData")

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
  mutate(Proportion_Trimmed_Kept = Non_Chimeras / sequence_counts_trimmed) %>%
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
  row.names = FALSE,
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

save.image(file = "../data/working/8_output.RData")
