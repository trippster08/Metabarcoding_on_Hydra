# /bin/bash

data=$(readlink -f data)
path_to_results=${data}/results
genes=()
truncation_values=()
possible_genes=(COI 18S MiFish 16S 28SAnth 16Sbac)
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

gene_num=${#genes[@]}
truncation_value_num={#truncation_values[@]}

if (( truncation_value_num < 2 * gene_num )); then
  echo "There are not enough truncation values given, there should be two truncation values given for each gene."
elif (( truncation_value_num > 2 * gene_num )); then
  echo "There are two many trunction values given, there should be two truncation values given for each gene."
else 
  echo "Here are the genes you are analyzing, with their respective truncation values:"
  for (( i=0; i<gene_num; i++ )); do
    echo "${genes[i]} R1: ${truncation_values[2*i]} R2: ${truncation_values[2*i+1]}"
  done
fi


for gene in ${genes}; do
  if [[ ! " ${possible_genes[@]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]]; then
    echo "Error: ${gene} is not a valid gene name. Valid genes are: ${possible_genes[@]}."
    exit 1
  fi
done
#echo ${gene_num}
#echo ${genes[@]}
#echo ${truncation_values[@]}

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