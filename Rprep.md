The first time you run this pipeline on Hydra, you need to install several R libraries first, and this needs to be done by hand. This may take some time, and it is recommended you use an interactive node.  

Start and interactive session. You can use more than one cpu to speed up installation of packages.
```
qrsh -pe mthread 4
```
Next, load the R module. There are two R modules currently installed on Hydra (4.4.0 and 4.4.1), with the default being 4.4.0. We want to use 4.4.1, so you have to specify that when you load the module. Installing these will take some time, especially tidyverse (>15 min) and dada2 (~ 10 min), so plan accordingly.
```
module load gcc/14.2.0
module load tools/R/4.4.1
module load tools/cmake

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
BiocManager::install("ShortRead", ask = FALSE, Ncpus = 4)
```

Install any other libraries you may need. Libraries will only need to be installed once. If you get a message saying some packages have more recent versions available,
and asking if you want to update them, chose "1: ALL".
There have been some issues with installing certain packages (such as DADA2, ShortRead, seqinr, etc). If you are getting errors when attempting to install these packages, see our [R and RStudio help page](https://confluence.si.edu/spaces/HPC/pages/385975502/Using+the+RStudio+Server). This page currently says it is only for the RStudio Server, but it also helps with command-line R problems.


```
install.packages("digest", Ncpus = 4)
install.packages("ggplot2", Ncpus = 4)
install.packages("dplyr", Ncpus = 4)
install.packages("tibble", Ncpus = 4)
install.packages("readr", Ncpus = 4)
install.packages("stringr", Ncpus = 4)
install.packages("tidyr", Ncpus = 4)
install.packages("seqinr", Ncpus = 4)
install.packages("ape", Ncpus = 4)
install.packages("R.utils", Ncpus = 4)
install.packages("data.table", Ncpus = 4)


```
After all packages have been restored, load them just to make sure they work.
```
library(dada2)
library(ShortRead)
library(digest)
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(tidyr)
library(seqinr)
library(ape)
library(R.utils)
library(data.table)
```
If all libraries load appropriately, quit R and continue with pipeline.
```
quit()
```