# Metabarcoding_on_Hydra

This protocol is for paired-end demultiplexed miseq sequences that have sufficient overlap to merge R1 and R2, and are going to be run on Hydra, the Smithsonian Institutions High Performance Cluster. This pipeline contains multiple .job files for submission to Hydra. Each job contains all the commands necessary to a point where results can be evaluated and a decision for parameters for the next job will needs to be made. For example, the first submission includes using Cutadapt to remove primers from raw reads and creating and visualizing quality plots from trimmed reads using DADA2 in R. At this point a decision must be

However, before analysing you data, you must make sure the necessary programs are installed, and the illumina demultiplexed sequences have been downloaded.
