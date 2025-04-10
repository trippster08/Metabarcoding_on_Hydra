# REMOVE CHIMERAS ##############################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## Remove Chimeric Sequences ===================================================

# Load the RData from "quality_plot_multigene.R"
load("../data/working/6_merge.RData")

# Here we remove chimera sequences. Use seqtabXXX if you removed long or short
# sequences above.
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

print(paste(
  "This is the number of ASVs for your chimera-free Sequence-Table:",
  length(colnames(seqtab_nochim)),
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
  sequences = as.list(repseq_chimera),
  names = repseq_chimera_md5,
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

# Save all the objects created to this point in this section
save.image(file = "../data/working/7_chimera.RData")

print("Job 7_chimera.job has finished")
