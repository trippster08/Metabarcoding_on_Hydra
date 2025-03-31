# TRIM FILTER DENOISE MERGE ####################################################
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
args <- commandArgs(trailingOnly = TRUE)
gene <- args[1]
path_to_trimmed <- args[2]
truncF <- as.numeric(args[3])
truncR <- as.numeric(args[4])

# Load the RData from "quality_plot_multigene.R"
load("../data/working/2_qual.RData")

# Next we are going to rename all our objects with the specific gene in the name
# (such as COI or 18S or 12S or whatever gene is assigned to "gene") to remove
# the gene name
object_names <- ls()
new_object_names <- gsub(paste0("_", gene), "", object_names)
for (i in seq_along(object_names)) {
  assign(new_object_names[i], get(object_names[i]))
}

## Trimmed Sequences ===========================================================
# This creates a vector of the path for forward reads (R1, called trimmed_F).
trimmed_F <- sort(list.files(
  path_to_trimmed,
  pattern = "_R1.fastq.gz",
  full.names = TRUE
))

trimmed_R <- sort(list.files(
  path_to_trimmed,
  pattern = "_R2.fastq.gz",
  full.names = TRUE
))

# Make a vector of sample names from your trimmed reads.
sample_names_trimmed <- sapply(
  strsplit(basename(trimmed_F), "_S\\d{1,3}_"),
  `[`,
  1
)
# Give the vectors names
names(trimmed_F) <- sample_names_trimmed
names(trimmed_R) <- sample_names_trimmed

## Get Read Counts of Trimmed Samples ==========================================
# Count the number of reads in each trimmed sample. Since cutadapt only
# keeps paired reads, we only need to count forward samples.
sequence_counts_trimmed <- sapply(trimmed_F, function(file) {
  fastq_data <- readFastq(file)
  length(fastq_data)
})
names(sequence_counts_trimmed) <- sample_names_trimmed

# This creates file paths for the reads that will be quality filtered with dada2
# in the next step.
filtered_F <- file.path(
  "../data/working",
  paste0(
    "filtered_sequences/",
    gene
  ),
  paste0(
    sample_names_trimmed,
    "_F_filt.fastq.gz"
  )
)
filtered_R <- file.path(
  "../data/working",
  paste0(
    "filtered_sequences/",
    gene
  ),
  paste0(
    sample_names_trimmed,
    "_R_filt.fastq.gz"
  )
)

# This inserts sample names to these newly created file paths.
names(filtered_F) <- sample_names_trimmed
names(filtered_R) <- sample_names_trimmed

# This filters all reads depending upon the quality (as assigned by the user)
# and trims the ends off the reads for all samples as determined by the quality
# plots. It also saves into "out" the number of reads in each sample that were
# filtered, and how many remain after filtering

# On Windows set multithread=FALSE. If you get errors while running this, change
# to multithread = FALSE, because "error messages and tracking are not handled
# gracefully when using the multithreading functionality".

# "truncLen=c(i,j)" is how you tell Dada2 where to truncate all forward (i) and
# reverse (j) reads. Using "0" means reads will not be truncated.
# maxEE sets how many expected errors are allowed before a read is filtered out.

# The amount to truncate is a common question, and very unsettled. I usually
# truncate at the point just shorter than where the red line (proportion of
# reads) in the quality plot reaches 100%.

# Most pipelines use a low maxEE (maximum number of expected errors), but I tend
# to relax this value (from 0,0 to 6,6) because it increases the number of reads
# that are kept, and Dada2 incorporates quality scores in its error models, so
# keeping poorer-quality reads does not adversely effect the results, except in
# very low quality reads. However, increasing maxEE does increase computational
# time.

filtered_summary <- filterAndTrim(
  trimmed_F,
  filtered_F,
  trimmed_R,
  filtered_R,
  truncLen = c(truncF, truncR),
  maxN = 0,
  maxEE = c(4, 4),
  rm.phix = TRUE,
  truncQ = 2,
  compress = TRUE,
  multithread = TRUE,
  verbose = TRUE
)

# After filtering, if you have any samples with no reads, you much remove them.
# This step changes filtered_F and filtered_R to only contain the names of
# samples with reads. Do this only if there are samples in
# "filtered_summary" with zero reads.
exists <- file.exists(filtered_F) & file.exists(filtered_R)
print(
  "Here are the number of samples that were removed because they no longer contain reads after filtering"
)
length(filtered_F) - length(filtered_F[exists])

filtered_F <- filtered_F[exists]
filtered_R <- filtered_R[exists]


# Set a path to the directory with the dada2-filtered reads.
path_to_filtered <- "../data/working/filtered_sequences"

# Get sample names for filtered reads
sample_names_filtered <- sapply(
  strsplit(basename(filtered_F), "_[FR]_filt"),
  `[`,
  1
)

# Count how many reads remain in each sample after filtering
sequence_counts_filtered <- sapply(filtered_F, function(file) {
  fastq_data <- readFastq(file)
  length(fastq_data)
})
# Name the counts with sample names
names(sequence_counts_filtered) <- sample_names_filtered

print(paste(
  "Here are the read counts for each filtered",
  gene,
  "sample:",
  sep = " "
))
sequence_counts_filtered

# Export filtered_summary as a tsv
write.table(
  filtered_summary,
  file = paste0(
    "../data/results/",
    gene,
    "/",
    project_name,
    "_filtered_reads_",
    gene,
    ".tsv"
  ),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = NA
)

# Save all the objects created to this point in this section
save.image(file = paste0("../data/working/", gene, "_3_filter.RData"))

print("Job 3_filter_multigene.job and this analysis has finished")
