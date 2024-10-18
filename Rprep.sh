# This loads the most recent version of R, currently 4.1.1.

module load tools/R/4.4.1

# Rscript is a way of running R commands (here through an .R file, but it can
# also run commands directly) without running R interactively
Rscript Rprep.R
