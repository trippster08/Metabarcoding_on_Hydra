# R AND HYDRA PREPARATION ######################################################

## Istall and Load R Libraries =================================================
# These are the libraries that will be used for this pipeline. Not all will be
# used for each persons particular project, but it does not hurt to have them
# all installed, loaded, and ready to be used when needed. You only need to
# install these one time.  When trying to install the first one, hydra will ask
# you if want to install into your personal library. Answer yes, and it will
# suggest a path to install all libraries into. Change the path if you want,
# or accept the default directory.
# Install  BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
# Install Dada2. You may get an error telling you to install a different version
# of Dada2. Change "3.19" to whatever version RStudio tells you.
BiocManager::install("dada2", version = "3.19")
# Install Phyloseq
BiocManager::install("phyloseq")
# Install DECIPHER
BiocManager::install("DECIPHER")

# Install other libraries you may need. Libraries will only need to be installed once.
install.packages("digest")
install.packages("tidyverse")
install.packages("seqinr")
install.packages("ape")
install.packages("filesstrings")
