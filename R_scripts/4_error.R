# LEARN ERRORS #################################################################
## Load Libraries ==============================================================
suppressMessages(library(dada2, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(digest, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(seqinr, warn.conflicts = FALSE, quietly = TRUE))
suppressMessages(library(ShortRead, warn.conflicts = FALSE, quietly = TRUE))

## File Housekeeping ===========================================================

# Load the RData from "quality_plot_multigene.R"
load("data/working/3_filter.RData")
## Error Estimation ============================================================

# Here we use a portion of the data to determine error rates. These error rates
# will be used in the next (denoising) step to narrow down the sequences to a
# reduced and corrected set of unique sequences

# Create a list to contain the gene-specific error rates
errors <- setNames(vector("list", length(genes)), genes)

# Make a loop to determine errors for each gene. First add read direction to
# each gene in the list, then model errors
for (gene in genes) {
  cat("\nModelling error rates and creating plots for", gene, ".\n")
  errors[[gene]] <- list(F = NULL, R = NULL)
  for (direction in c("F", "R")) {
    errors[[gene]][[direction]] <- learnErrors(
      filtered_reads[[gene]][[direction]],
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
}

# We can visualize the estimated error rates to make sure they don't look too
# crazy. The red lines are error rates expected under the "...nominal defintion
# of the Q-score." The black dots are "...observed error rates for each
# consensus quality score." The black line shows the "...estimated error rates
# after convergence of the machine-learning algorithm." I think the main things
# to look at here are to make sure that each black line is a good fit to the
# observed error rates, and that estimated error rates decrease with increased
# quality.

# First we make a list to store the gene-specific error plots
error_plots <- setNames(vector("list", length(genes)), genes)
# Make a loop to create error plots for each gene. First add read direction to
# each gene in the list, then make plots, then export plots as a pdf
for (gene in genes) {
  error_plots[[gene]] <- list(F = NULL, R = NULL)

  for (direction in c("F", "R")) {
    err <- errors[[gene]][[direction]]
    error_plots[[gene]][[direction]] <- plotErrors(err, nominalQ = TRUE)
    ggsave(
      filename = file.path(
        path_to_working,
        paste0(
          project_name,
          "_",
          gene,
          "_errorplot_",
          direction,
          ".pdf"
        )
      ),
      plot = error_plots[[gene]][[direction]],
      width = 9,
      height = 9
    )
  }
}

# Save all the objects created to this point in this section
save.image(file = "data/working/4_error.RData")

cat("\n4_error.job is complete, error rates evaluated, and plots created.\n")
