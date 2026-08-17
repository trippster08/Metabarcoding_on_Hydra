#!/bin/bash
path_to_data=$(realpath ./data)
project_name=$(basename "$PWD")
gene="$@"

qsub -o logs/tax_rdp_${gene}_hydra.log -N ${project_name}_${gene}_tax_rdp \
jobs/9_assign_taxonomy_dada2_loop.job ${gene}