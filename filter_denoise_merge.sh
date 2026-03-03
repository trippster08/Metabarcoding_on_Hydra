#!/bin/bash

path_to_data=$(realpath ./data)
path_to_results=${path_to_data}/results
genes=()
truncation_values=()
available_genes=()

# Loop through all arguments
for arg in "$@"; do
  if [[ $arg =~ ^[0-9]+$ ]]; then
    # If the argument is a number, add it to the numbers array
    truncation_values+=("$arg")
  else
    # If the argument is not a number, add it to the gene names array
    genes+=("$arg")
  fi
done
# assign the number of variables in the genes array as gene_num (the number of genes)
gene_num=${#genes[@]}
# assign the number of variables in the truncation_values array as truncation_value_num
truncation_value_num=${#truncation_values[@]}

# Loop through all the files in the primer folder and get the names from
# each forward primer file. Then add them to the variable of available primer
# names
for file in primers/*F.fas; do
  primer_file=${file##*/}
  gene_name=${primer_file%F.fas}
  available_genes+=(${gene_name})
done
#echo ${available_genes[@]}
#echo ${genes[@]}
# Check to make sure there are two truncation values for every gene value
if (( truncation_value_num < 2 * gene_num )); then
  echo "There are not enough truncation values given, there should be two truncation values given for each gene."
  exit 1
elif (( truncation_value_num > 2 * gene_num )); then
  echo "There are two many trunction values given, there should be two truncation values given for each gene."
  exit 1
else # If there are, print out the genes and their respective truncation values
  echo "Here are the genes you are analyzing, with their respective truncation values:"
  for (( i=0; i<gene_num; i++ )); do
    echo "${genes[i]} R1: ${truncation_values[2*i]} R2: ${truncation_values[2*i+1]}"
  done
fi

# Check to make sure the gene names given are correct
for gene in ${genes[@]}; do
  # If one of the variables is not a valid primer, print error and list of primers
  if [[ ! " ${available_genes[@]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]]; then
    echo "Error: ${gene} is not a valid gene name. Valid genes are: ${available_genes[@]}."
    exit 1
  fi
done
#echo ${gene_num}
#echo ${genes[@]}
#echo ${truncation_values[@]}



# For each gene, submit its own job(s)
for (( i=0; i<gene_num; i++ )); do
  gene=${genes[i]}
  r1=${truncation_values[2*i]}
  r2=${truncation_values[2*i+1]}

  echo "Processing gene: ${gene} (R1 trunc=${r1}, R2 trunc=${r2})"


  # Gene-specific checkpoint paths (you'll need R to write these)
  filter_rdata="data/working/${gene}_3_filter.RData"
  error_rdata="data/working/${gene}_4_error.RData"
  denoise_rdata="data/working/${gene}_5_denoise.RData"
  merge_rdata="data/working/${gene}_6_merge.RData"
  chimera_rdata="data/working/${gene}_7_chimera.RData"
  output_rdata="data/working/${gene}_8_output.RData"


  if [[ ! -f "$filter_daata" ]]; then
    # Step 3: filter, per gene
    qsub -o "logs/filter_${gene}.log" -N "filter_${gene}" \
      jobs/3_filter.job "$gene" "$r1" "$r2"
    echo "Submitted filter job for ${gene}"

  elif [[ ! -f "$error_rdata" ]]; then
    # Step 4: error modelling, per gene
    qsub -o "logs/error_${gene}.log" -N "error_${gene}" \
      jobs/4_error.job "$gene"
    echo "Filtering complete for ${gene}; starting error modelling."

  elif [[ ! -f "$denoise_rdata" ]]; then
    # Step 5: denoise, per gene
    qsub -o "logs/denoise_${gene}.log" -N "denoise_${gene}" \
      jobs/5_denoise.job "$gene"
    echo "Error modelling complete for ${gene}; starting denoising."

  elif [[ ! -f "$merge_rdata" ]]; then
    # Step 6: merge, per gene
    qsub -o "logs/merge_${gene}.log" -N "merge_${gene}" \
      jobs/6_merge.job "$gene"
    echo "Denoising complete for ${gene}; starting merging."

  elif [[ ! -f "$chimera_rdata" ]]; then
    # Step 7: chimera
    qsub -o "logs/chimera_${gene}.log" -N "chimera_${gene}" \
      jobs/7_chimera.job "$gene"
    echo "Merging complete for ${gene}; starting chimera removal."

  elif [[ ! -f "$output_rdata" ]]; then
    # Step 8: output
    qsub -o "logs/output_${gene}.log" -N "output_${gene}" \
      jobs/8_output.job "$gene"
    echo "Chimera removal complete for ${gene}; exporting results."

  else
    echo "All steps for gene ${gene} have already completed."
  fi

done