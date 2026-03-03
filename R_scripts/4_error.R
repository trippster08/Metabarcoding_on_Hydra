# LEARN ERRORS #################################################################
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
    stop(paste0("No gene argument provided in job file 4_error_", gene, ".job."))
}

# Save gene name and truncation values as three separate variables
gene <- args[1]
# Load the RData from "3_filter.RData"
load(paste0("data/working/3_filter_", gene, ".RData"))
## Error Estimation ============================================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences

cat("\nModelling error rates and creating plots for", gene, "\n")
# Create list for forward and reverse error rates
errors <- list(F = NULL, R = NULL)
for (direction in c("F", "R")) {
  errors[[direction]] <- learnErrors(
    filtered_reads[[direction]],
    nbases = 1e+08,
    errorEstimationFunction = loessErrfun,
    multithread = TRUE,
    randomize = FALSE,
    MAX_CONSIST = 10,
    OMEGA_C = 0,
    qualityType = "Auto",
    verbose = FALSE
  )
}

# We can visualize the estimated error rates to make sure they don't look too
# crazy. The red lines are error rates expected under the "...nominal defintion
# of the Q-score." The black dots are "...observed error rates for each
# consensus quality score." The black line shows the "...estimated error rates
# after convergence of the machine-learning algorithm." I think the main things
# to look at here are to make sure that each black line is a good fit to the
# observed error rates, and that estimated error rates decrease with increased
# quality.

# Create list for forward and reverse error plots
error_plots <- list(F = NULL, R = NULL)
# Make a loop to create error plots.
for (direction in c("F", "R")) {
  err <- errors[[direction]]
  error_plots[[direction]] <- plotErrors(err, nominalQ = TRUE)
  ggsave(
    filename = file.path(
      path_to_results,
      paste0(
        "plots/",
        project_name,
        "_",
        gene,
        "_errorplot_",
        direction,
        ".pdf"
      )
    ),
    plot = error_plots[[direction]],
    width = 9,
    height = 9
  )
}


# Save all the objects created to this point in this section
save.image(file = paste0("data/working/4_error_", gene, ".RData"))

cat("\n4_error.job for", gene, "is complete, error rates evaluated, and plots created.\n")
