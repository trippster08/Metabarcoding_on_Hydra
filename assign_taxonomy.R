suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(rBLAST, warn.conflicts = FALSE, quietly = TRUE))


## Trim Reads ==================================================================
numcores <- Sys.getenv("NSLOTS")
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

setwd(paste0("/scratch/genomics/macdonaldk/edna_workshop_2025/", sample))
load("data/results/feattab.RData")


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
