# R AND HYDRA PREPARATION ######################################################

## Istall and Load R Libraries =================================================
# These are the libraries that will be used for this pipeline. Not all will be
# used for each persons particular project, but it does not hurt to have them
# all installed, loaded, and ready to be used when needed.
# Install  BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
# Install Dada2. You may get an error telling you to install a different version
# of Dada2. Change "3.17" to whatever version RStudio tells you.
BiocManager::install("dada2", version = "3.17")
# Install Phyloseq
BiocManager::install("phyloseq")
# Install DECIPHER
BiocManager::install("DECIPHER")

# Install and other libraries you may need (or install through
# "Install Packages" window). Libraries will only need to be installed once.
install.packages("digest")
install.packages("tidyverse")
install.packages("seqinr")
install.packages("ape")
install.packages("filesstrings")
