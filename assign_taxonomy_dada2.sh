#!/bin/bash
path_to_data=$(realpath ./data)
project_name=$(basename "$PWD")
genes="$@"

for gene in ${genes}; do
  qsub -o logs/tax_dada2_${gene}_hydra.log -N ${project_name}_${gene}_tax_dada2 \
  jobs/9_assign_taxonomy_dada2_loop.job ${gene} ${path_to_data} ${project_name}
  sleep 0.1
done
