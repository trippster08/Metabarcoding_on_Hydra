The first time you run this pipeline on Hydra, you need to install several R libraries first, and this needs to be done by hand. This may take some time, and it is recommended you use an interactive node.  

Start and interactive session. You can use more than one cpu to speed up installation of packages.
```
qrsh -pe mthread 4
```
Next, load the R module. There are two R modules currently installed on Hydra (4.4.0 and 4.4.1), with the default being 4.4.0. We want to use 4.4.1, so you have to specify that when you load the module. Installing these will take some time, especially tidyverse (>15 min) and dada2 (~ 10 min), so plan accordingly.
```
module load tools/R/4.4.1
```
Next we need to start R. To do that, you just type "R".  You can also go to an RStudio server and open this page and run the installations through the server.

```
R
```
Install BiocManager. If this is the first R library you have install, you will get an error message saying that you don't have permission, and asking if you want to install it in your personal library. Respond with `yes`.
```
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
```
Install all of the libraries needed through BiocManager.
```
BiocManager::install("dada2", ask = FALSE, Ncpus = 4)
BiocManager::install("phyloseq", ask = FALSE, Ncpus = 4)
BiocManager::install("msa", ask = FALSE, Ncpus = 4)
BiocManager::install("DECIPHER", ask = FALSE, Ncpus = 4)
BiocManager::install("ShortRead", ask = FALSE, Ncpus = 4)
BiocManager::install("rBLAST", ask = FALSE, Ncpus = 4)
```

Install any other libraries you may need. Libraries will only need to be installed once. If you get a message saying some packages have more recent versions available,
and asking if you want to update them, chose "1: ALL".
```
install.packages("digest", Ncpus = 4)
install.packages("tidyverse", Ncpus = 4)
install.packages("seqinr", Ncpus = 4)
install.packages("ape", Ncpus = 4)
install.packages("vegan", Ncpus = 4)
install.packages("patchwork", Ncpus = 4)
install.packages("remotes", Ncpus = 4)
install.packages("R.utils", Ncpus = 4)
install.packages("phylotools", Ncpus = 4)
install.packages("data.table", Ncpus = 4)

remotes::install_github("ropensci/taxize", upgrade = TRUE)
remotes::install_github("fkeck/refdb", upgrade = TRUE)
remotes::install_github("tobiasgf/lulu", upgrade = TRUE)
remotes::install_github("boldsystems-central/BOLDconnectR", upgrade = TRUE)
install.packages("rMSA", repos = "https://mhahsler.r-universe.dev")
```
After all packages have been restored, load them just to make sure they work.
```
library(dada2)
library(tidyverse)
library(digest)
library(seqinr)
library(ape)
library(DECIPHER)
library(phyloseq)
library(msa)
library(ShortRead)
library(rBLAST)
library(vegan)
library(patchwork)
library(phylotools)
library(data.table)
```
If all libraries load appropriately, quit R and continue with pipeline.
```
q()
```
