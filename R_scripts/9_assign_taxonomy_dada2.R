# ASSIGN TAXONOMY ##############################################################
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
  stop("No gene argument provided in job file assign_taxonomy_dada2_loop.job")
}

# get gene name from argument
gene <- args[1]


load(paste0(
  "data/working/8_output_",
  gene,
  ".RData"
))


# You need a reference database. This one is a COI one that works for most
# marine COI sequences, but has all insects removed. You can provide your own
# but it needs to be in the correct DADA2 format (see
# https://benjjneb.github.io/dada2/training.html) for formatting or downloading
# DADA2 formatted databases.

# If you are using your own reference database, put the path to the database
# here. This path is the most recent silva database for SSR bacteria and
#archaea. You can also use your own database, but it needs to be in the correct
# DADA2 format (see https://benjjneb.github.io/dada2/training.html) for
# formatting or downloading 
reference_fasta = "/scratch/nmnh_lab/macdonaldk/ref/silva_nr99_v138.2_toSpecies_trainset.fa.gz"

# Make sure the taxonomic levels match the levels in your reference database.
# These are the levels for the silva database.
tax_levels <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "species"
)

# Assign taxonomy
taxonomy <- dada2::assignTaxonomy(
  seqtab_nochim,
  reference_fasta,
  taxLevels = tax_levels,
  tryRC = FALSE,
  minBoot = 50,
  outputBootstraps = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

# Convert ASV to ASV-md5 hash
taxonomy_md5 <- tibble::as_tibble(taxonomy$tax, rownames = "ASV") %>%
  dplyr::mutate(md5 = repseq_nochim_md5_asv$md5)
boot_md5 <- tibble::as_tibble(taxonomy$boot, rownames = "ASV") %>%
  dplyr::mutate(md5 = repseq_nochim_md5_asv$md5)



# Save this object, it took a long time to get.
save(
  taxonomy,
  file = paste0(
    "data/working/9_taxonomy_rdp_",
    gene,
    ".Rdata"
  )
)

write.table(
  taxonomy_md5,
  file = paste0(
    "data/results/",
    gene,
    "/",
    project_name,
    "_dada2_taxonomy_tax.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write.table(
  boot_md5,
  file = paste0(
    "data/results/",
    gene,
    "/",
    project_name,
    "_dada2_taxonomy_boot.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
