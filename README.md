![GitHub Release](https://img.shields.io/github/v/release/trippster08/Metabarcoding_on_Hydra?color=green)
[![DOI](https://zenodo.org/badge/849911002.svg)](https://doi.org/10.5281/zenodo.15635236)

# Metabarcoding_on_Hydra

This protocol is for paired-end demultiplexed miseq sequences that have sufficient overlap to merge R1 and R2, and are going to be run on Hydra, the Smithsonian Institution's High Performance Cluster. This pipeline contains multiple .job files for submission to Hydra. These jobs are called by two shell scripts. The first shell script submits a job to remove primers from raw reads and create and visualize quality plots from trimmed reads, at which point trimming parameters need to be decided. The second shell script submits a job to filter, denoise, and remove chimeras, and outputs multiple files, including a Feature-Table, Sequence-Table, and a representative sequence fasta. This pipeline is also designed to secondarily demuliplex by primer if mulitple genes were sequences in a single run.

This pipeline is designed run on [Hydra](https://confluence.si.edu/display/HPC/High+Performance+Computing), Smithsonian's HPC, using [cutadapt](https://github.com/marcelm/cutadapt/) to remove primers from raw Illumina reads, and [DADA2](https://benjjneb.github.io/dada2/) in R to filter, trim, denoise, remove chimeras, and merge trimmed reads.   The pipeline assumes you have a current hydra account and are capable of accessing the SI network, either in person or through VPN. Our pipeline is specifically written for MacOS, but is compatible with Windows. See [Hydra on Windows PCs](https://confluence.si.edu/display/HPC/Logging+into+Hydra) for differences between MacOS and Windows in accessing Hydra.

___
**NOTICE!**

The directions have changed since the last version of this pipeline. Please **carefully** read the directions below.
___

## Hydra Configuration 
All programs will be run through pre-installed Hydra modules/programs, so there is no need for the user to change any Hydra configurations or install any programs.

### Log into Hydra
Open the terminal app and log onto Hydra. You will need your hydra account password.
```
ssh USERNAME@hydra-login01.si.edu
```
 or
```
ssh USERNAME@hydra-login02.si.edu
```
### Project-specific Directory 
Go to the the directory assigned to you for short-term storage of large data-sets. Typically this will be `/scratch/genomics/USERNAME/`. Replace USERNAME with your hydra username.
```
cd /scratch/genomics/USERNAME
```
Make a project-specific directory, with the subdirectory `data/raw/`. -p allows you to create subdirectories and any parental ones that don't already exist (in this case, PROJECT). Again, replace PROJECT with your project name. Everything will be run from this main project directory.

```
mkdir -p PROJECT/data/raw

```
### Transfer Files to Hydra 
From the main project directory, download the pipeline using `wget` (see code below). This downloads a compressed file that contains all job files (\*.job), shell scripts (\*.sh), R scripts (\*.R), and primer definition files (\*.fas) necessary for your analysis. 

```
wget https://github.com/trippster08/Metabarcoding_on_Hydra/archive/refs/heads/main.zip
```
Copy and paste the script below into your terminal. It will unzip the pipeline, and move all the \*.sh files and the job, R_script and primer directories from your newly unzipped directory into the main project directory. It will also delete the now-empty pipeline directory and zipped download.
```
unzip main.zip
mv Metabarcoding_on_Hydra-main/*.sh \
Metabarcoding_on_Hydra-main/jobs \
Metabarcoding_on_Hydra-main/primers \
Metabarcoding_on_Hydra-main/R_scripts .
rm -r Metabarcoding_on_Hydra-main main.zip

```
### Get Raw Reads
Your raw reads should be copied into `data/raw/`. Download to your local computer and use scp or filezilla to upload to `data/raw/`. See [Transferring Files to/from Hydra](https://confluence.si.edu/pages/viewpage.action?pageId=163152227) for help with transferring files between Hydra and your computer. *NOTE: Remove any "undetermined" read files from the folder containing your raw reads. You do not want to include these reads in your analyses.*
If your files are saved as a dropbox link to a .zip file, you can upload your raw reads directly to hydra without first saving to your computer (remembering that /scratch is not for long-term storage of data) using wget. However, for large datasets we often run into errors both in uploading the file and in unzipping it. If you want to attempt this, I can walk you through the process.

### Preparing R
The first time you run this pipeline, you may need to install libraries that may be needed: (`DADA2`, `tidyverse`, `seqinr`, `ape`, `digest`, etc). Currently, you have to do this manually, please see [Installing R libraries on Hydra](https://github.com/trippster08/Metabarcoding_on_Hydra/blob/main/Rprep.md) for directions on how to install these libraries.

## Running the Pipeline
The basic pipeline runs as two shell scripts, each of which submits one or more jobs to Hydra and results in a Feature-table of ASVs and associated files. The first shell script submits a job to trim primer regions off all reads using Cutadapt, and demultiplex reads by primer pair if mixed samples are run on one machine. This first job also creates quality plots for trimmed reads using DADA2 in R. If trimming has already occured, the first shell script will skip trimming and submit only the job to create quality plots.  After examining the quality plots, you chose truncation values for the DADA2 filtering step and run the second shell script. This script submits a job that runs DADA2 in R, using your chosen truncation values to filter trimmed reads, then denoises, merges, and removes chimeras before outputing a feature-table, representative-sequence fasta, and various other results. If you have demultiplexed by gene in the trimming step, the job will run for all genes and save your results in a gene-specific directory. If your job fails or gets interupted at any point, just restart this second script and it will determine the last completed step and restart the pipeline at the beginning of the next step. Scripts for data visualization, taxonomic assignment and other outputs are or will be added as independent modules with their own directions.

You can monitor job progress by examining the logs in the logs/ directory. If there are no errors, there should only be two logs for a completed analysis: cutadapt.log and filter.log. If errors occurred and the jobs needed to be restarted, additional logs may be created.

### Trim Primers and Create Quality Plots
This pipeline uses Cutadapt to remove primers from reads, remove poly-G tails if necessary, and remove short reads. It will also demultiplex by amplicon primer if more than one gene region was part of the same run. The following shell script submits a job to run Cutadapt for all samples in a run on Hydra. Run this script, followed by the names of the one or more gene regions in the run. As an example, `sh trim_and_quality_plot.sh COI` will search for and remove "mlCOIintF/jgCOIR" from each pair of reads, while `sh trim_and_quality_plot.sh MiFish COI` will search for both the COI primers and MiFish 12S prmers. It will also demultiplex reads, moving each read into a gene-specific folder depending on the primers removed. In our example the gene-specific folders are `data/working/trimmed_sequences/COI` and `data/working/trimmed_sequences/MiFish`. Following trimming, gene-specific quality plots are created using DADA2

We currently have primer sets for five gene regions, COI, 18S-v4, MiFish-12S, 16S-v4 for soil microbiomes, and 28S for Anthozoans. The default names for these five regions are "COI", "18S", "MiFish", "16Sbac" and "28SAnth. If you give it an incorrect or unavailable name, you will get a message telling you so and giving a list of compatible names. If you have additional primers to trim, I can add them to the default primer list, or you can add custom primers yourself on a per-analysis basis.

#### Custom Primers
You can currently trim custom primers in addition to those that come with the pipeline. See the instructions at the bottom of [Custom Primers](https://github.com/trippster08/Metabarcoding_on_Hydra/blob/main/primers/primer_info.md) for directions on how to add your own primer files with sequences. You then append the name of your custom primers (i.e. gene region) to `sh trim_and_quality_plot.sh`.


```
sh trim_and_quality_plot.sh <gene_region_used>
```
or, for more than one gene
```
sh trim_and_quality_plot.sh <gene1_region_used> <gene2_region_used> <gene2_region_used> 
```
#### Examine Quality Plots
Depending on how you are accessing Hydra, either look at your quality plots in your viewer or download and open locally. For these plots, the green line is the mean quality score at that position, the orange lines are the quartiles (solid for median, dashed for 25% and 75%) and the red line represents the proportion of reads existing at that position. Determine where you want to truncate reads based on quality and read number. The amount to truncate is a common question, and very unsettled. I usually truncate at the point just shorter than where the red line (proportion of reads) in the quality plot drops precipitously from 100%. Never trim beyond this point, you will lose most of your samples.

### Filter, Trim and Merge Using DADA2
After trimming and examining the quality plots, we run a second shell script to submit a job to run DADA2 in R. This job will filter, quality trim, and truncate reads, followed by denoising, merging, and removing chimeric ASVs. The output will be multiple files placed in `data/results/`. Results include a Feature-Table (rows of ASV's as md5 hashs and columns of samples), a Sample-Table (rows of samples and columns of ASV's as md5 hashs), a rep-seq fasta file (ASV's with their respective md5 hash as a name), a two-column table of ASVs and md5 hashs, and a table that tracks read counts through the process. You submit this job through a shell script followed by the truncation values for R1 and R2. For example, `sh filter_denoise_merge.sh COI 260 260` will input trimmed COI reads and truncate all reads at 260 bp in the filtering and trimming step. If you have multiple genes, input can be either "gene1 gene2 R1_gene1 R2_gene2 R1_gene2 R2_gene2" or "gene1 R1_gene1 R2_gene1 gene2 R1_gene2 R2_gene2". The main necessity is that truncation values are in the same order as genes and R1 comes before R2.
```
sh filter_denoise_merge.sh <gene_region_used> <R1_truncation_value> <R2_truncation_value>
```
or
```
sh filter_denoise_merge.sh <gene1_region_used> <R1_gene1_truncation_value> <R2_gene1_truncation_value> <gene2_region_used> <R1_gene2_truncation_value> <R2_gene2_truncation_value>
```
This completes the basic pipeline for metabarcoding illumina runs, resulting in gene-specfic feature-tables and represenatative-sequence fastas. As mentioned earlier, taxonomic assignment and further data visualization will be added as seperate modules soon (hopefully).