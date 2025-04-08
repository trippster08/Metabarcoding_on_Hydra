# ASSIGN TAXONOMY ##############################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(rBLAST, warn.conflicts = FALSE, quietly = TRUE))


## File Housekeeping ===========================================================
args <- commandArgs(trailingOnly = TRUE)
ref <- args[1]


# Load the RData from "quality_plot_multigene.R"
load("../data/working/8_output.RData")

taxonomy <- assignTaxonomy(
  seqtab_nochim,
  "/scratch/nmnh_lab/macdonaldk/ref/midori_COI_genus_dada2.fasta",
  taxLevels = c(
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "species"
  ),
  tryRC = FALSE,
  minBoot = 50,
  outputBootstraps = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

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
