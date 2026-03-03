# REMOVE CHIMERAS ##############################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Get argument containing the gene name from job file.
args <- commandArgs(trailingOnly = TRUE)

# Check to make sure there is an argument for the gene name
if (length(args) < 1) {
    stop(paste0("No gene argument provided in job file 7_chimera_", gene, ".job"))
}

# get gene name from argument
gene <- args[1]

# Load the RData from "6_merge_<gene>.RData"
load(paste0("data/working/6_merge_", gene, ".RData"))

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences.
cat(
  "\nRemoving chimeric sequences and creating new sequence table for",
  gene,
  "\n"
)
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# Print to the log the number of ASVs remaining in the table.
cat(
  "\nThis is the number of ASVs for your chimera-free",
  gene,
  "Sequence-Table:",
  length(colnames(seqtab_nochim)),
  "\n"
)

# Make a list of the ASVs that are considered chimeras, in case you want to
# look at them later
chimeras_list <- isBimeraDenovoTable(
  seqtab,
  multithread = TRUE,
  verbose = TRUE
)

# This makes a new vector containing all the ASV's (unique sequences) returned
# by dada2.
repseq_all <- getSequences(seqtab)
# Get a list of just chimera ASVs by filtering all sequences by chimera list
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

# Export chimeric sequences as fastas
write.fasta(
  sequences = as.list(repseq_chimera),
  names = repseq_chimera_md5,
  open = "w",
  as.string = FALSE,
  file.out = file.path(
    path_to_results,
    paste0(
      "additional_results/",
      project_name,
      "_",
      gene,
      "_rep-seq_chimeras.fas"
    )
  )
)
## Examine Sequence Lengths ==================================================

# This shows the length of the representative sequences (ASV's). Typically,
# there are a lot of much longer and much shorter sequences.

# Count the number of bp for each gene-specific ASV from the sequence-tables
seq_length_table <- table(nchar(getSequences(seqtab_nochim)))

# Export this table as a .tsv
write.table(
  seq_length_table,
  file = file.path(
    path_to_results,
    paste0(
      "additional_results/",
      project_name,
      "_ASV_lengths_table_",
      gene,
      ".tsv"
    )
  ),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

# Save all the objects created to this point in this section
save.image(file = paste0("data/working/7_chimera_", gene, ".RData"))

cat("\n7_chimera.job is done and chimeric sequences removed for gene:", gene, "\n")
