The first time you run this pipeline on Hydra, you need to install several R libraries first, and this needs to be done by hand. First, load the R module. There are two R modules currently installed on Hydra (4.4.0 and 4.4.1), with the default being 4.4.0. We want to use 4.4.1, so you have to specify that when you load the module. Installing these will take some time, especially tidyverse (>15 min) and dada2 (~ 10 min), so plan accordingly.
```
module load tools/R/4.4.1
```
Next we need to start R. To do that, you just type "R".  

```
R
```
Install BiocManager. If this is the first R library you have install, you will get an error message saying that you don't have permission, and asking if you want to install it in your personal library. Respond with `yes`.
```
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
```
Install DADA2. You may get an error telling you to install a different version. Change "3.19" to whatever version it tells you.
```
BiocManager::install("dada2", version = "3.19")
```

Install the rest of the libraries. You will be asked to choose a mirror for some of the libraries. Pick whichever you prefer, I dont think it matters.
```
install.packages("digest")
install.packages("tidyverse")
install.packages("seqinr")
install.packages("ape")
```
After all packages have been restored, load them just to make sure they work.
```
library(dada2)
library(tidyverse)
library(digest)
library(seqinr)
library(ape)
```
If all libraries load appropriately, quit R and continue with pipeline.
```
q()
```
