#!/bin/bash
path_to_data=$(realpath ./data)
project_name=$(basename "$PWD")
genes="$@"
for gene in ${genes}; do
  qsub -o logs/${gene}_blast_hydra.log -N ${project_name}_${gene}_blast \
  jobs/assign_taxonomy_blast_loop.job ${gene} ${path_to_data} ${project_name}
  sleep 0.1
done