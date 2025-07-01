# LEARN ERRORS #################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================
args <- commandArgs(trailingOnly = TRUE)
gene_num <- as.numeric(args[1])
genes <- args[2:(gene_num + 1)]

# Load the RData from "quality_plot_multigene.R"
load("data/working/3_filter.RData")
## Error Estimation ============================================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences
errors <- setNames(vector("list", length(genes)), genes)

for (gene in genes) {
  errors[[gene]] <- list(F = NULL, R = NULL)

  errors[[gene]]$F <- learnErrors(
    filtered_reads[[gene]]$F,
    nbases = 1e+08,
    errorEstimationFunction = loessErrfun,
    multithread = TRUE,
    randomize = FALSE,
    MAX_CONSIST = 10,
    OMEGA_C = 0,
    qualityType = "Auto",
    verbose = FALSE
  )
  errors[[gene]]$R <- learnErrors(
    filtered_reads[[gene]]$R,
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
error_plots <- setNames(vector("list", length(genes)), genes)
for (gene in genes) {
  error_plots[[gene]] <- list(F = NULL, R = NULL)

  for (direction in c("F", "R")) {
    err = errors[[gene]][[direction]]
    error_plots[[gene]][[direction]] <- plotErrors(err, nominalQ = TRUE)
    ggsave(
      filename = file.path(
        path_to_results,
        paste0(
        project_name,
        "_errorplots",
        direction,
        "_",
        gene,
        ".pdf"
      ),
      plot = error_plots[[gene]][[direction]],
      width = 9,
      height = 9
      )
    )
  }
}


# Save all the objects created to this point in this section
save.image(file = "data/working/4_error.RData")

print("Job 4_error_multigene.job and this analysis has finished")
