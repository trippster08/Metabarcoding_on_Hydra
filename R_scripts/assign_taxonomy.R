# ASSIGN TAXONOMY ##############################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(rBLAST, warn.conflicts = FALSE, quietly = TRUE))


## File Housekeeping ===========================================================

# You need a reference database. This one is a COI one that works for most
# marine COI sequences, but has all insects removed. You can provide your own
# but it needs to be in the correct DADA2 format (see
# https://benjjneb.github.io/dada2/training.html) for formatting or downloading
# DADA2 formatted databases.

# If you are using your own reference database, put the path to the database
# here
reference_fasta = "/scratch/nmnh_lab/macdonaldk/metabarcoding/ref/MIDORI2_UNIQ_NUC_GB260_CO1_DADA2_noInsect_unzipped.fasta"

# Load the RData from "quality_plot_multigene.R"
load("../data/working/8_output.RData")
# Make sure the taxonomic levels match the levels in your reference database
tax_levels <- c(
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "species"
)

# Assign taxonomy
taxonomy <- assignTaxonomy(
  seqtab_nochim,
  reference_fasta,
  taxLevels = tax_levels,
  tryRC = FALSE,
  minBoot = 50,
  outputBootstraps = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

# Save this object, it took a long time to get.
save(taxonomy, file = "data/working/tax_rdp.Rdata")

write.table(
  taxonomy$tax,
  file = paste0("data/results/", sample, "taxonomy_tax.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write.table(
  taxonomy$boot,
  file = paste0("data/results/", sample, "taxonomy_boot.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
