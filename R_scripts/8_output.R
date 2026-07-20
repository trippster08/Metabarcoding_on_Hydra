# OUTPUT RESULTS ###############################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tibble, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Get argument containing the gene name from job file.
args <- commandArgs(trailingOnly = TRUE)

# Check to make sure there is an argument for the gene name
if (length(args) < 1) {
  stop("No gene argument provided in job file 8_output_<gene>.job")
}

# get gene name from argument
gene <- args[1]
# Load the RData from "7_chimera_<gene>.RData" 
load(paste0("data/working/7_chimera_", gene, ".RData"))

## Track Reads Through Dada2 Process ===========================================

# Here, we look at how many reads made it through each step. This is similar to
# the stats table that we look at in Qiime2 or "track" from the DADA2 tutorial.
# This is a good quick way to see if something is wrong (i.e. only a small
# proportion make it through).

# First make a table for the post-filtered samples, including denoised,
# merged, and non-chimeric read counts

getN <- function(x) sum(dada2::getUniques(x))
sequence_counts_postfiltered <- tibble::tibble(
  Sample_ID = sample_names_filtered,
  Denoised_Reads_F = sapply(denoised$F, getN),
  Denoised_Reads_R = sapply(denoised$R, getN),
  Merged_Reads = sapply(merged_reads, getN),
  Non_Chimeras = as.integer(rowSums(seqtab_nochim))
)

  # Then we are going to add the postfiltered read count data to the three count
  # data objects we already have (raw, trimmed, filtered).
track_reads <- tibble::tibble(
  Sample_ID = names(sequence_counts_raw),
  Raw_Reads = as.numeric(sequence_counts_raw),
) %>%
  dplyr::left_join(
    tibble::tibble(
      Sample_ID = names(sequence_counts_trimmed[[gene]]),
      Trimmed_Reads = as.numeric(sequence_counts_trimmed[[gene]])
    ),
    join_by(Sample_ID)
  ) %>%
  dplyr::left_join(
    tibble::tibble(
      Sample_ID = names(sequence_counts_filtered),
      Filtered_Reads = as.numeric(sequence_counts_filtered)
    ),
    join_by(Sample_ID)
  ) %>%
  dplyr::left_join(
    sequence_counts_postfiltered,
    join_by(Sample_ID)
  ) %>%
  dplyr::mutate(Proportion_Trimmed_Passed = Non_Chimeras / Trimmed_Reads) %>%
  dplyr::mutate(Proportion_Raw_Passed = Trimmed_Reads / Raw_Reads) %>%
  dplyr::select(
    Sample_ID,
    Raw_Reads,
    Trimmed_Reads,
    Filtered_Reads,
    Denoised_Reads_F,
    Denoised_Reads_R,
    Merged_Reads,
    Non_Chimeras,
    Proportion_Trimmed_Passed,
    Proportion_Raw_Passed
  )

# Export this table as a .tsv
  write.table(
    track_reads,
    file = file.path(
      path_to_results,
      paste0(
        "additional_results/",
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

## Export Sequence-Table =======================================================
# This exports a sequence-table: columns of ASV's, rows of samples, and
# values = number of reads.

# If you have mulitple Miseqruns for the same project that will need to be
# combined for further analyses, you may want to name this file
# "PROJECTNAME_MISEQRUN_sequence-table.tsv" to differentiate different runs.
# In "5 Metabarcoding_R_Pipeline_RStudio_ImportCombine" we'll show how to
# combine data from separate runs for analyses.

# NOTE!!!
# The Sequence-Table we have now (seqtab_nochim) is very unwieldy, since each
# column name is an entire ASV. Instead, we will convert ASVs using the md5
# encryption model, creating a 32bit representative "hash" of each ASV. Every
# hash is essentially unique to the ASV it is representing. We would then
# replace the ASVs in the column headings with their representative md5 hash.
# However, having an ASV hash as a column heading requires the creation of a
# Representative Sequence list, which tells us which hash represents which ASV.

### Create And Use md5 Hash ----------------------------------------------------
# To create a Sequence list with md5 hash instead of ASVs, we first need to
# create a list of md5 hash's of all ASV's.

# Use the program digest (in a For Loop) to create a new vector containing the
# unique md5 hashes of the representative sequences (ASV's). This results in
# identical feature names to those assigned in Qiime2.

repseq_nochim <- dada2::getSequences(seqtab_nochim)

repseq_nochim_md5 <- c()
for (i in seq_along(repseq_nochim)) {
  repseq_nochim_md5[i] <- digest::digest(
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
  file = file.path(
    path_to_results,
    paste0(
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
repseq_nochim_md5_asv <- tibble::tibble(
  md5 = repseq_nochim_md5,
  ASV = repseq_nochim
)

## Create and Export Feature-Table ===========================================
# This creates and exports a feature-table: row of ASV's (shown as a md5 hash
# instead of sequence), columns of samples, and values = number of reads. With
# this table you will also need a file that relates each ASV to it's
# representative md5 hash. We download this in the next section.

# Transpose the sequence-table, and convert the result into a tibble.
feattab_nochim_md5 <- t(seqtab_nochim_md5)

write.table(
  feattab_nochim_md5,
  file = file.path(
    path_to_results,
    paste0(
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

## Export Representative Sequences table/fasta ===============================
# Here we export our our representative sequences, either as a fasta (with the
# md5 hash as the ASV name), or as a table with ASV and md5 hash as columns.

# This exports all the ASVs in fasta format, with ASV hash as the sequence
# name. This is analogous to the representative sequence output in Qiime2.
seqinr::write.fasta(
  sequences = as.list(repseq_nochim_md5_asv$ASV),
  names = repseq_nochim_md5_asv$md5,
  open = "w",
  as.string = FALSE,
  file.out = file.path(
    path_to_results,
    paste0(
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
  repseq_nochim_md5_asv,
  file = file.path(
    path_to_results,
    paste0(
      project_name,
      "_representative-seq-md5-table_",
      gene,
      ".tsv"
    )
  ),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

## Create and Export Quality Scores =========================================
# Create a table to get mean sample quality (phred) scores, and % of reads
# within each sample with phred scores > Q30

### Create a Quality Score function -------------------------------------------
# o look at each fastq.gz, get the sum of scores and total
# number of bases for that .fastq.gz
fastq_stats <- function(file) {
  # Create a connection to the .fastq.gz file. "rt" means read text. This does
  # not unzip the file. We would just use file instead of gzfile for a .fastq
  con <- gzfile(file, "rt")
  # This will close the file after the function is done, even if the function
  # creashes
  on.exit(close(con))
  # Totals scores and bases...so you divide the total score by the number of
  # bases, and we have the average score for the file
  total_score <- 0
  total_bases <- 0
  q30_bases <- 0
  # Create a loop that goes on until we stop it (infinite loop)
  repeat {
    # read the lines, 400000 at a time, which equals 100000 reads (each read
    # contains 4 lines). This reduces memory usage 
    lines <- readLines(con, n = 400000)
    # when we run out of lines, rreadLines = 0, and we stop the loop
    if (length(lines) == 0)
      break
    # Take each 4th line, which is the quality score, make a vector of them
    qual <- lines[seq(4, length(lines), 4)]
    # save phred score from ASCI quality. This takes each symbol, converts
    # to its ASCI numerical value, and subtracts 33, and does it for each read
    qvals <- unlist(
      lapply(
        qual,
        function(x) as.integer(charToRaw(x)) - 33L
      ),
      use.names = FALSE
    )

    # total all phred scores
    total_score <- total_score + sum(qvals)
    # count number of bases
    total_bases <- total_bases + length(qvals)
    # count number reads with pred scores greater than 30
    q30_bases <- q30_bases + sum(qvals >= 30)
  }
  
  c(
    # calculate mean phred and percentage of reads with phred over 30
    mean_phred = total_score / total_bases,
    pct_Q30 = 100 * q30_bases / total_bases
  )
}

### Create and Export Quality Score Table ------------------------------------
all_reads <- c(filtered_reads$F, filtered_reads$R)
direction <- c(
  rep("F", length(filtered_reads$F)),
  rep("R", length(filtered_reads$R))
)

stopifnot(is.character(all_reads))
stopifnot(all(file.exists(all_reads)))

results <- t(sapply(all_reads, fastq_stats))

quality_summary <- data.frame(
  file = basename(all_reads),
  mean_phred = results[, "mean_phred"],
  pct_Q30 = results[, "pct_Q30"],
  row.names = NULL
)

quality_summary <- quality_summary[
  order(quality_summary$file),
]

write.table(
  quality_summary,
  file = file.path(
    path_to_results,
    paste0(
      "/additional_results/",
      project_name,
      "_quality_scores_",
      gene,
      ".tsv"
    )
  ),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

## Create and Export feature-to-fasta ========================================
# This creates a fasta file containing all the ASV's for each sample. Each ASV
# will be labeled with the sample name, ASV hash, and number of reads of that
# ASV in that sample. This was derived from a python script from Matt Kweskin
# called featuretofasta.py (hence the name).

### Create a Sequence-List Table -----------------------------------------------
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

seqtab_nochim_tall <- tibble::as_tibble(
  seqtab_nochim,
  rownames = "sample"
) %>%
  tidyr::pivot_longer(
    !sample,
    names_to = "ASV",
    values_to = "count"
  ) %>%
  subset(count != 0)

# Save the ASV sequences from the sequence-list table
# (seqtab_nochim_tall_nozera) as a new list.
repseq_tall <- seqtab_nochim_tall$ASV

# Convert the sequences into md5 hashs, as we did earlier. md5 hashs are
# consistent across jobs, meaning identical sequences from different projects
# or being converted by different programs will result in the same hash (i.e.
# hashs here will match hashs above)
repseq_tall_md5 <- c()
for (i in seq_along(repseq_tall)) {
  repseq_tall_md5[i] <- digest::digest(
    repseq_tall[i],
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
seqtab_nochim_tall_md5 <- seqtab_nochim_tall %>%
  dplyr::mutate(md5 = repseq_tall_md5) %>%
  dplyr::mutate(sample_md5_count = paste(sample, md5, count, sep = "_")) %>%
  dplyr::select(sample, md5, count, sample_md5_count, ASV)

### Export feature-to-fasta fas file -----------------------------------------
# Create a fasta-formatted file of each row sequence (i.e. ASV), with a
# heading of "sample_feature_count".
seqinr::write.fasta(
  sequences = as.list(seqtab_nochim_tall_md5$ASV),
  names = seqtab_nochim_tall_md5$sample_md5_count,
  open = "w",
  as.string = FALSE,
  file.out = file.path(
    path_to_results,
    paste0(
      "additional_results/",
      project_name,
      "_feature-to-fasta_",
      gene,
      ".fas"
    )
  )
)

save.image(file = paste0("data/working/8_output_", gene, ".RData"))

cat("\n8_output.job and this analysis for", gene, "has finished.\n")
