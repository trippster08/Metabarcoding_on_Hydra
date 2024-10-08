# Metabarcoding_on_Hydra

This protocol is for paired-end demultiplexed miseq sequences that have sufficient overlap to merge R1 and R2, and are going to be run on Hydra, the Smithsonian Institutions High Performance Cluster. This pipeline contains multiple .job files for submission to Hydra. Each job contains all the commands necessary to a point where results can be evaluated and a decision for parameters for the next job needs to be made. For example, the first submission includes using Cutadapt to remove primers from raw reads and creating and visualizing quality plots from trimmed reads using DADA2 in R, at which point trimming parameters need to be decided.

This pipeline is designed run multiple samples simultanteously on [Hydra](https://confluence.si.edu/display/HPC/High+Performance+Computing), Smithsonian's HPC, using [cutadapt](https://github.com/marcelm/cutadapt/) to remove primer from raw Illumina reads, and [DADA2](https://benjjneb.github.io/dada2/) in R to filter, trim, and denoise trimmed reads.   The pipeline assumes you have a current hydra account and are capable of accessing the SI network, either in person or through VPN. Our pipeline is specifically written for MacOS, but is compatible with Windows. See [Hydra on Windows PCs](https://confluence.si.edu/display/HPC/Logging+into+Hydra) for differences between MacOS and Windows in accessing Hydra.

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
Transfer the pipeline to `jobs/` in your Hydra account. This downloads a compressed file that contains all job files (`*.job`), shell scripts (`*.sh`), R scripts (`*.R`), and primer definition files (`*.fas`) necessary for your analysis. This command downloads a compressed file that will become a directory upon unzipping.

```
cd jobs
wget https://github.com/trippster08/Metabarcoding_on_Hydra/archive/refs/heads/main.zip
```
Unzip the pipeline, and move all the shell, .job, and .R files from your newly unzipped directory into the job directory and the primer folder into the main project directory. Delete the now-empty pipeline directory.
```
unzip main.zip
mv Metabarcoding_on_Hydra-main/* .
mv primers ..
rm -r Metabarcoding_on_Hydra-main
```

Your raw reads should be copied into `data/raw/`. We are working on a way to load raw read files directly from their dropbox folder into our Hydra `data/raw/` directory. In the meantime, download to your local computer and use scp or filezilla to upload to `data/raw/`. See https://confluence.si.edu/pages/viewpage.action?pageId=163152227 for help with transferring files between Hydra and your computer. 

## Running the Pipeline



