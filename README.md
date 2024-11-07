# Metabarcoding_on_Hydra

This protocol is for paired-end demultiplexed miseq sequences that have sufficient overlap to merge R1 and R2, and are going to be run on Hydra, the Smithsonian Institutions High Performance Cluster. This pipeline contains multiple .job files for submission to Hydra. Each job contains all the commands necessary to a point where results can be evaluated and a decision for parameters for the next job needs to be made. For example, the first submission includes using Cutadapt to remove primers from raw reads and creating and visualizing quality plots from trimmed reads using DADA2 in R, at which point trimming parameters need to be decided.

This pipeline is designed run on [Hydra](https://confluence.si.edu/display/HPC/High+Performance+Computing), Smithsonian's HPC, using [cutadapt](https://github.com/marcelm/cutadapt/) to remove primer from raw Illumina reads, and [DADA2](https://benjjneb.github.io/dada2/) in R to filter, trim, and denoise trimmed reads.   The pipeline assumes you have a current hydra account and are capable of accessing the SI network, either in person or through VPN. Our pipeline is specifically written for MacOS, but is compatible with Windows. See [Hydra on Windows PCs](https://confluence.si.edu/display/HPC/Logging+into+Hydra) for differences between MacOS and Windows in accessing Hydra.

## Local Computer Configuration 
Make a project directory, and mulitple subdirectories on your local computer. Make this wherever you want to store your projects. Hydra is not made for long-term storage, so raw sequences, jobs, results, etc should all be kept here when your analyses are finished. Although it is not necessary, I use the same directory pattern locally as I use in Hydra. 

Make sure to replace "PROJECT" with your project name throughout.
```
mkdir -p PROJECT/data/raw PROJECT/jobs
```
Your raw reads should be in `data/raw/`. 
NOTE: As currently designed, this pipeline has a few naming requirements for raw reads. Reads should be `fastq.gz` or `fastq` formated, and needs to start with a unique sample name (that contains no underscores) and read number (either R1 or R2) later in the filename. Both sample name and read number must be separated from the remainder of the filename with underscores. Also, Hydra does not allow jobs names to start with a number, so if your sample names start with a number, change the name by adding at least one letter to the beginning of the filename (I usually use the initials of the researcher) before running this pipeline.

## Hydra Configuration 
All programs will be run through pre-installed Hydra modules, so there is no need for the user to change any Hydra configurations or install any programs.

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
Make a project-specific directory, with the following subdirectories: `jobs/` and `data/raw/`. -p allows you to create subdirectories and any parental ones that don't already exist (in this case, PROJECT). I use the same directory tree here as on my local computer, to lessen confusion. Again, replace PROJECT with your project name.
This pipeline is not dependent upon the directory tree shown, so you can set up your space differently, if you prefer. The only two directories that are required are `/data` and `jobs/` but you can name them whatever you like, and neither necessarily have to be in any particular place.This pipeline does create seveal new directories, including `data/working/`, `data/results/`, and `jobs/logs/`. If you don't want these directories created, or want them in different places, they can be changed in the shell scripts. 

```
mkdir -p PROJECT/data/raw PROJECT/jobs

```
### Transfer Files to Hydra 
Transfer the pipeline to `jobs/` in your Hydra account. This downloads a compressed file that contains all job files (`*.job`), shell scripts (`*.sh`), R scripts (`*.R`), and primer definition files (`*.fas`) necessary for your analysis. This command downloads a compressed file that will become a directory upon unzipping. Don't forget to move into your jobs folder first: `cd PROJECT/jobs`.

```
wget https://github.com/trippster08/Metabarcoding_on_Hydra/archive/refs/heads/main.zip
```
Unzip the pipeline, and move all the shell, .job, and .R files from your newly unzipped directory into the job directory and the primer folder into the main project directory. Delete the now-empty pipeline directory.
```
unzip main.zip
mv Metabarcoding_on_Hydra-main/* .
mv primers ..
rm -r Metabarcoding_on_Hydra-main
```

Your raw reads should be copied into `data/raw/`. We are working on a way to load raw read files directly from their dropbox folder into our Hydra `data/raw/` directory. In the meantime, download to your local computer and use scp or filezilla to upload to `data/raw/`. See [Transferring Files to/from Hydra](https://confluence.si.edu/pages/viewpage.action?pageId=163152227) for help with transferring files between Hydra and your computer. 

### Preparing R
The first time you run this pipeline, you may need to install libraries that may be needed: (`DADA2`, `tidyverse`, `seqinr`, `ape`, and `digest`). Currently, you have to do this manually, please see [Installing R libraries on Hydra](https://github.com/trippster08/Metabarcoding_on_Hydra/blob/main/Rprep.md) for directions on how to install these libraries.

## Running the Pipeline
This pipeline runs as two shell scripts, each of which submits a job to Hydra. The first shell script submits a job to trim primer regions off all reads, and demultiplex reads by primer pair if mixed samples are run on one machine, using Cutadapt. This first job also creates quality plots for trimmed reads using DADA2 in R.  After examining the quality plots, you chose truncation values for the DADA2 filtering step and run the second shell script. This script submits a job that runs DADA2 in R, using your chosen truncation values to filter trimmed reads, then denoises, merges, and removes chimeras before outputing a feature-table and representative-sequence fasta.

### Trim Primers and Create Quality Plots for Single-Gene Runs
This pipeline uses Cutadapt to remove primers from reads, remove poly-G tails if necessary, and remove short reads. The following shell script submits a job to run Cutadapt for all samples in a run on Hydra. Run this script, followed by the path to the raw reads and the gene region used. As an example, `sh trim_and_quality_plot.sh /scratch/genomics/USER/PROJECT/data/raw COI` will search for and remove "mlCOIintF/jgCOIR" from each pair of reads.
Following trimming, quality plots are created using DADA2

We currently have primer sets for three gene regions, COI, 18S-v4, and MiFish-12S. The default names for these three regions are "COI", "18S", and "12S", but the script can take some variations of each gene name, and if you use one that is not compatible, it will give you a list of compatible names.

```
sh trim_and_quality_plot.sh <path_to_raw_reads> <gene_region_used>
```
### Trim Primers and Create Quality Plots for Multiple-Gene Runs
You can currently also analyse runs that contain reads two gene regions. Cutadapt will separate reads based on the primer-pair removed, and move each trimmed read into a gene-specific trimmed-reads directory. To run on multi-gene runs, run the following shell script, followed by the path to the raw reads and both gene regions used. As an example, `sh trim_and_quality_plot.sh /scratch/genomics/USER/PROJECT/data/raw COI 18S` will search for and remove "mlCOIintF/jgCOIR" and "18S_V4_F/18S_V4_R" from each pair of reads, and move those reads into trimmed_reads/COI/ or trimmed_reads/18S/ depending upon which primer pair was removed. Quality plots will then be created for these demultiplexed Gene-specific reads and saved in data/results/.

We currently have primer sets for three gene regions, COI, 18S-v4, and MiFish-12S. The default names for these three regions are "COI", "18S", and "12S", but the script can take some variations of each gene name, and if you use one that is not compatible, it will give you a list of compatible names.

```
sh trim_and_quality_plot_multigene.sh <path_to_raw_reads> <gene1_region_used> <gene2_region_used>
```
#### Examine Quality Plots
Depending on your how you are accessing Hydra, either look at your quality plots in your viewer or download. # For these plots, the green line is the mean quality score at that position, the orange lines are the quartiles (solid for median, dashed for 25% and 75%) and the red line represents the proportion of reads existing at that position. Determine where you want to truncate reads based on quality and read number. The amount to truncate is a common question, and very unsettled. I usually truncate at the point just shorter than where the red line (proportion of reads) in the quality plot drops precipitously from 100%.

### Filter, Trim and Merge Using DADA2 for Single-Gene runs
After trimming and examining the quality plots, we submit a second job to run DADA2 in R. This job will filter, quality trim, and truncate reads, followed by denoising, merging, and removing chimeric ASVs. This will output three final results files, placing them in data/results/: a Feature-Table (rows of ASV's as md5 hashs and columns of samples), a rep-seq fasta file (ASV's with their respective md5 hash as a name), and a two-column table of ASVs and md5 hashs. You submit this job through a shell script followed by the directory of the trimmed reads and truncation values for R1 and R2. For example, `sh filter_denoise_merge.sh /scratch/genomics/USER/PROJECT/data/working/trimmed_reads/ 260 260` will input trimmed reads and truncate all reads at 260 bp in the filtering and trimming step.
```
sh filter_denoise_merge.sh <path_to_trimmed_reads> <R1_truncation_value> <R2_truncation_value>
```
### Filter, Trim and Merge Using DADA2 for Multiple-Gene runs
After trimming and examining the quality plots, we submit a second job to run DADA2 in R. This job will filter, quality trim, and truncate reads, followed by denoising, merging, and removing chimeric ASVs. It will do this seperately for each gene region. This will output three final results files for each gene, and place them in a gene-specific results directory: a Feature-Table (rows of ASV's as md5 hashs and columns of samples), a rep-seq fasta file (ASV's with their respective md5 hash as a name), and a two-column table of ASVs and md5 hashs. Currently, the submission method is clunky and is a work in progress. You submit this job through a shell script followed by the directories of each gene region (gene1 and gene2), then the truncation values for gene1 R1 and R2, then the truncation values for gene2 R1 and R2. For example, `sh filter_denoise_merge.sh /scratch/genomics/USER/PROJECT/data/working/trimmed_reads/COI sh filter_denoise_merge.sh /scratch/genomics/USER/PROJECT/data/working/trimmed_reads/18S 260 260 250 250` will input trimmed reads from both gene-specific directory and truncate all reads for COI at 260 bp and for 18S at 250 bp.
```
sh filter_denoise_merge.sh <path_to_trimmed_reads_gene1> <path_to_trimmed_reads_gene2> <R1_truncation_value_gene1> <R2_truncation_value_gene1> <R1_truncation_value_gene2> <R2_truncation_value_gene2>
```
