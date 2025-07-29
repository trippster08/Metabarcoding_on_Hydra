# /bin/bash

path_to_data=$(cd ./data && pwd)
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
# ecah forward primer file. Then add them to the variable of available primer
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

# Check to see if the .RData file for each job has been created. If it has, then that
# job has finished, and submit the next job. The first job for which the .RData file
# does not exist, start that job.
if [[ ! -f "data/working/3_filter.RData" ]]; then
  qsub -o logs/filter.log -N filter \
    jobs/3_filter.job ${gene_num} ${genes[@]} ${truncation_values[@]}
elif [[ ! -f "data/working/4_error.RData" ]]; then
  qsub -o logs/error.log -N error \
    jobs/4_error.job ${gene_num} ${genes[@]}
  echo "Filtering is already complete; we will start with error modelling."
elif [[ ! -f "data/working/5_denoise.RData" ]]; then
  qsub -o logs/denoise.log -N denoise \
    jobs/5_denoise.job ${gene_num} ${genes[@]}
  echo "Error modelling is already complete; we will start with denoising."
elif [[ ! -f "data/working/6_merge.RData" ]]; then
  qsub -o logs/merge.log -N merge \
    jobs/6_merge.job ${gene_num} ${genes[@]}
  echo "Denoising is already complete; we will start with merging reads."
elif [[ ! -f "data/working/7_chimera.RData" ]]; then
  qsub -o logs/chimera.log -N chimera \
    jobs/7_chimera.job ${gene_num} ${genes[@]}
  echo "Merging of reads is already complete; we will start with chimera removal."
elif [[ ! -f "data/working/8_output.RData" ]]; then
  qsub -o logs/output.log -N output \
    jobs/8_output.job ${gene_num} ${genes[@]}
  echo "Chimera removal is already complete; we will start with exporting results."
else
  echo "All steps for this analysis have already completed"
fi