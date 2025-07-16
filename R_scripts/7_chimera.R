# REMOVE CHIMERAS ##############################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
# Load the RData from "quality_plot_multigene.R"
load("data/working/6_merge.RData")

## Remove Chimeric Sequences ===================================================

# Here we remove chimera sequences.

# Make three lists that will be populated later by gene-specific data
seqtab_nochim <- setNames(vector("list", length(genes)), genes)
repseq_chimera_md5 <- setNames(vector("list", length(genes)), genes)
seq_length_table <- setNames(vector("list", length(genes)), genes)
# Loop through each gene, remove chimeras from the gene-specific sequence-table
for (gene in genes) {
  seqtab_nochim[[gene]] <- removeBimeraDenovo(
    seqtab[[gene]],
    method = "consensus",
    multithread = TRUE,
    verbose = TRUE
  )

  # Print to the log the number of ASVs remaining in the table.
  print(paste(
    "This is the number of ASVs for your chimera-free",
    gene,
    "Sequence-Table:",
    length(colnames(seqtab_nochim[[gene]])),
    sep = " "
  ))

  # Make a list of the ASVs that are considered chimeras, in case you want to
  # look at them later
  chimeras_list <- isBimeraDenovoTable(
    seqtab[[gene]],
    multithread = TRUE,
    verbose = TRUE
  )

  # This makes a new vector containing all the ASV's (unique sequences) returned
  # by dada2.
  repseq_all <- getSequences(seqtab[[gene]])
  # Get a list of just chimera ASVs by filtering all sequences by chimera list
  repseq_chimera <- repseq_all[chimeras_list]
  # Make and add md5 hash to the repseq_chimera
  repseq_chimera_md5[[gene]] <- c()
  for (i in seq_along(repseq_chimera)) {
    repseq_chimera_md5[[gene]][i] <- digest(
      repseq_chimera[i],
      serialize = FALSE,
      algo = "md5"
    )
  }

  # Export chimeric sequences as fastas
  write.fasta(
    sequences = as.list(repseq_chimera),
    names = repseq_chimera_md5[[gene]],
    open = "w",
    as.string = FALSE,
    file.out = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
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
  seq_length_table[[gene]] <- table(nchar(getSequences(seqtab_nochim[[gene]])))

  # Export this table as a .tsv
  write.table(
    seq_length_table[[gene]],
    file = file.path(
      path_to_results[[gene]],
      paste0(
        "/",
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
}
# Save all the objects created to this point in this section
save.image(file = "data/working/7_chimera.RData")

print("Job 7_chimera.job is done and chimeric sequences removed")
