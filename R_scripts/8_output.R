# OUTPUT RESULTS ###############################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================

# Load the RData from "quality_plot_multigene.R"
load("data/working/7_chimera.RData")

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2 or "track" from the DADA2 tutorial.
# This is a good quick way to see if something is wrong (i.e. only a small
# proportion make it through).

# First make a table for the post-filtered samples, including denoised,
# merged, and non-chimeric read counts
sequence_counts_postfiltered <- setNames(vector("list", length(genes)), genes)
track_reads <- setNames(vector("list", length(genes)), genes)
for (gene in genes) {
  getN <- function(x) sum(getUniques(x))
  sequence_counts_postfiltered[[gene]] <- tibble(
    Sample_ID = sample_names_filtered[[gene]],
    Denoised_Reads_F = sapply(denoised[[gene]]$F, getN),
    Denoised_Reads_R = sapply(denoised[[gene]]$R, getN),
    Merged_Reads = sapply(merged_reads[[gene]], getN),
    Non_Chimeras = as.integer(rowSums(seqtab_nochim[[gene]]))
  )

  # Then we are going to add the postfiltered read count data to the three count
  # data objects we already have (raw, trimmed, filtered).
  track_reads[[gene]] <- tibble(
    Sample_ID = names(sequence_counts_raw),
    Raw_Reads = as.numeric(sequence_counts_raw),
  ) %>%
    left_join(
      tibble(
        sequence_counts_trimmed[[gene]],
        Sample_ID = names(sequence_counts_trimmed[[gene]]),
        Trimmed_Reads = as.numeric(sequence_counts_trimmed[[gene]])
      ),
      join_by(Sample_ID)
    ) %>%
    left_join(
      tibble(
        sequence_counts_filtered[[gene]],
        Sample_ID = names(sequence_counts_filtered[[gene]]),
        Filtered_Reads = as.numeric(sequence_counts_filtered[[gene]])
      ),
      join_by(Sample_ID)
    ) %>%
    left_join(
      sequence_counts_postfiltered[[gene]],
      join_by(Sample_ID)
    ) %>%
    mutate(Proportion_Trimmed_Kept = Non_Chimeras / Trimmed_Reads) %>%
    mutate(Proportion_Gene = Trimmed_Reads / Raw_Reads) %>%
    select(
      Sample_ID,
      Raw_Reads,
      Trimmed_Reads,
      Filtered_Reads,
      Denoised_Reads_F,
      Denoised_Reads_R,
      Merged_Reads,
      Non_Chimeras,
      Proportion_Trimmed_Kept,
      Proportion_Gene
    )
}
# Export this table as a .tsv
for (gene in genes) {
  write.table(
    track_reads[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_track_reads_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )
}

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
for (gene in genes) {
  write.table(
    seqtab_nochim[[gene]],
    file = file.path(
      path_to_working,
      paste0(
        project_name,
        "_sequence-table_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )
}
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
#  table you will later use for your analyses (e.g. seqtab.nochim)
repseq_nochim <- setNames(vector("list", length(genes)), genes)
repseq_nochim_md5 <- setNames(vector("list", length(genes)), genes)
seqtab_nochim_md5 <- setNames(vector("list", length(genes)), genes)
repseq_nochim_md5_asv <- setNames(vector("list", length(genes)), genes)
feattab_nochim_md5 <- setNames(vector("list", length(genes)), genes)
# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.
for (gene in genes) {
  repseq_nochim[[gene]] <- getSequences(seqtab_nochim[[gene]])

  repseq_nochim_md5[[gene]] <- c()
  for (i in seq_along(repseq_nochim[[gene]])) {
    repseq_nochim_md5[[gene]][i] <- digest(
      repseq_nochim[[gene]][i],
      serialize = FALSE,
      algo = "md5"
    )
  }

  # Add md5 hash to the sequence-table from the DADA2 analysis.

  seqtab_nochim_md5[[gene]] <- seqtab_nochim[[gene]]
  colnames(seqtab_nochim_md5[[gene]]) <- repseq_nochim_md5[[gene]]

  # Export this sequence table with column headings as md5 hashs instead of ASV
  # sequences
  write.table(
    seqtab_nochim_md5[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_sequence-table-md5_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE,
    col.names = NA
  )

  # Create an md5/ASV table, with each row as an ASV and it's representative md5
  # hash.
  repseq_nochim_md5_asv[[gene]] <- tibble(
    md5 = repseq_nochim_md5[[gene]],
    ASV = repseq_nochim[[gene]]
  )

  ## Create and Export Feature-Table =============================================
  # This creates and exports a feature-table: row of ASV's (shown as a md5 hash
  # instead of sequence), columns of samples, and values = number of reads. With
  # this table you will also need a file that relates each ASV to it's
  # representative md5 hash. We download this in the next section.

  # Transpose the sequence-table, and convert the result into a tibble.
  feattab_nochim_md5[[gene]] <- t(seqtab_nochim_md5[[gene]])

  write.table(
    feattab_nochim_md5[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_feature-table_md5_",
        gene,
        ".tsv"
      )
    ),
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
    sequences = as.list(repseq_nochim_md5_asv[[gene]]$ASV),
    names = repseq_nochim_md5_asv[[gene]]$md5,
    open = "w",
    as.string = FALSE,
    file.out = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_rep-seq_",
        gene,
        ".fas"
      )
    )
  )

  # This exports all the ASVs and their respective md5 hashes as a two-column
  # table.
  write.table(
    repseq_nochim_md5_asv[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_representative_sequence_md5_table_",
        gene,
        ".tsv"
      )
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}
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
# The sequence-table has a column with sample names, and N columns of ASV's
# containing count values. We want all the count data to be in a single column,
# so we use a tidyr command called "pivot_longer" to make the table "tall",
# which means the table goes from 17 x 2811 to 47770 x 3 for example
# (47770 = 2810 x 17. 2810 instead of 2811 because the first column of the
# original table contains sample names, not counts). This makes the table tidier
# (meaning that each column is now a true variable).
seqtab_nochim_tall <- setNames(vector("list", length(genes)), genes)
repseq_tall <- setNames(vector("list", length(genes)), genes)
repseq_tall_md5 <- setNames(vector("list", length(genes)), genes)
seqtab_nochim_tall_md5 <- setNames(vector("list", length(genes)), genes)

for (gene in genes) {
  seqtab_nochim_tall[[gene]] <- as_tibble(
    seqtab_nochim[[gene]],
    rownames = "sample"
  ) %>%
    pivot_longer(
      !sample,
      names_to = "ASV",
      values_to = "count"
    ) %>%
    subset(count != 0)

  ## Create and Export feature-to-fasta ==========================================
  # This creates a fasta file containing all the ASV's for each sample. Each ASV
  # will be labeled with the sample name, ASV hash, and number of reads of that
  # ASV in that sample. This was derived from a python script from Matt Kweskin
  # called featuretofasta.py (hence the name).

  # Save the ASV sequences from the sequence-list table
  # (seqtab_nochim_tall_nozera) as a new list.
  repseq_tall[[gene]] <- seqtab_nochim_tall[[gene]]$ASV

  # Convert the sequences into md5 hashs, as we did earlier. md5 hashs are
  # consistent across jobs, meaning identical sequences from different projects or
  # being converted by different programs will result in the same hash (i.e.
  # hashs here will match hashs above)
  repseq_tall_md5[[gene]] <- c()
  for (i in seq_along(repseq_tall[[gene]])) {
    repseq_tall_md5[[gene]][i] <- digest(
      repseq_tall[[gene]][i],
      serialize = FALSE,
      algo = "md5"
    )
  }

  # Attach the ASV hashes as a column (called "md5") to the tall table. The
  # table should now have 4 columns, and each row of the "md5" column should
  # be a md5 hash of its respective ASV.
  # Also create a new column in this table that contains "sample", "feature",
  # and "count", concatenated. This is the heading for each sequence in the
  # fasta file created by Matt Kweskin's script "featuretofasta.py"
  seqtab_nochim_tall_md5[[gene]] <- seqtab_nochim_tall[[gene]] %>%
    mutate(md5 = repseq_tall_md5[[gene]]) %>%
    mutate(sample_md5_count = paste(sample, md5, count, sep = "_")) %>%
    select(sample, md5, count, sample_md5_count, ASV)

  # Create a fasta-formatted file of each row sequence (i.e. ASV), with a
  # heading of "sample_feature_count".
  write.fasta(
    sequences = as.list(seqtab_nochim_tall_md5[[gene]]$ASV),
    names = seqtab_nochim_tall_md5[[gene]]$sample_md5_count,
    open = "w",
    as.string = FALSE,
    file.out = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
        project_name,
        "_feature-to-fasta_",
        gene,
        ".fas"
      )
    )
  )
}
save.image(file = paste0("data/working/8_output.RData"))

print("Job 8_output.job and this analysis has finished")
