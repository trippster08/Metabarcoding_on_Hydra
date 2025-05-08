 # /bin/bash
trimmed=$(readlink -f ../data/working/trimmed_sequences)
gene1="$1"
gene2="$2"
truncF1="$3"
truncR1="$4"
truncF2="$5"
truncR2="$6"
COI=(COI coi CO1 co1 cox1 COX1)
MiFish=(12S MiFish mifish Mifish 12S_mifish 12S_MiFish 12s)
V4=(18S l8s V4 v4)

echo "Your trimmed reads are in this directory: ${trimmed}"
# Function to check if a term is in a group
is_in_group() {
  local term="$1"
  shift
  local group=("$@")
  for item in "${group[@]}"; do
    if [[ "$term" == "$item" ]]; then
      return 0
    fi
  done
  return 1
}

# Check if gene1 and gene2 are in any of the groups
if (is_in_group "$gene1" "${COI[@]}" || is_in_group "$gene1" "${MiFish[@]}" || is_in_group "$gene1" "${V4[@]}") &&
   (is_in_group "$gene2" "${COI[@]}" || is_in_group "$gene2" "${MiFish[@]}" || is_in_group "$gene2" "${V4[@]}"); then
  echo "You are denoising " ${gene1} "and " ${gene2}
else
  echo "gene1 and/or gene2 are not valid terms. Please use COI, MiFish, or V4 to define your two genes."
fi

# Check if there are read files in both the trimmed sequences directories
if
  [[ -z "$(ls ${trimmed}/${gene1}/*.fastq.gz 2>/dev/null | grep fastq)" ]] || \
  [[ -z "$(ls ${trimmed}/${gene2}/*.fastq.gz 2>/dev/null | grep fastq)" ]]  
then  
  echo "At least one of the directories holding your trimmed reads is empty"
  exit
fi

if
  [[ ! ${truncF1} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncR1} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncF2} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncR2} =~ ^[0-9]+$ ]]
then
  echo "Length to truncation values for one of your genes is not entered correctly. It must be a number between 0 (no trimming) and the length of your reads. Please enter <gene1> <gene2> <R1_truncation_value_for_gene1> <R2_truncation_value_for_gene1> <R1_truncation_value_for_gene2> <R2_truncation_value_for_gene2>'."
  exit
fi

#echo ${gene1}
#echo ${gene2}
#echo ${trimmed}
#echo ${truncF1}
#echo ${truncR1}
#echo ${truncF2}
#echo ${truncR2}

if [[ ! -f "../data/working/${gene1}_3_filter.RData" ]]; then
  qsub -o logs/filter_${gene1}.log -N filter_${gene1} \
    3_filter_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
elif [[ ! -f "../data/working/${gene1}_4_error.RData" ]]; then
  qsub -o logs/error_${gene1}.log -N error_${gene1} \
    4_error_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
  echo ${gene1} "filtering is already complete; we will start with error modelling."
elif [[ ! -f "../data/working/${gene1}_5_denoise.RData" ]]; then
  qsub -o logs/denoise_${gene1}.log -N denoise_${gene1} \
    5_denoise_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
  echo ${gene1} "error modelling is already complete; we will start with denoising."
elif [[ ! -f "../data/working/${gene1}_6_merge.RData" ]]; then
  qsub -o logs/merge_${gene1}.log -N merge_${gene1} \
    6_merge_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
  echo ${gene1} "denoising is already complete; we will start with merging reads."
elif [[ ! -f "../data/working/${gene1}_7_chimera.RData" ]]; then
  qsub -o logs/chimera_${gene1}.log -N chimera_${gene1} \
    7_chimera_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
  echo ${gene1} "merging is already complete; we will start with chimera removal."
elif [[ ! -f "../data/working/${gene1}_8_output.RData" ]]; then
  qsub -o logs/output_${gene1}.log -N output_${gene1} \
    8_output_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
  echo ${gene1} "chimera removal is already complete; we will start with exporting results."
else
  echo "All steps for ${gene1} have already completed"
fi

if [[ ! -f "../data/working/${gene2}_3_filter.RData" ]]; then
  qsub -o logs/filter_${gene2}.log -N filter_${gene2} \
    3_filter_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF1} ${truncR1}
elif [[ ! -f "../data/working/${gene2}_4_error.RData" ]]; then
  qsub -o logs/error_${gene2}.log -N error_${gene2} \
    4_error_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF1} ${truncR1}
  echo ${gene2} "filtering is already complete; we will start with error modelling."
elif [[ ! -f "../data/working/${gene2}_5_denoise.RData" ]]; then
  qsub -o logs/denoise_${gene2}.log -N denoise_${gene2} \
    5_denoise_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF1} ${truncR1}
  echo ${gene2} "error modelling is already complete; we will start with denoising."
elif [[ ! -f "../data/working/${gene2}_6_merge.RData" ]]; then
  qsub -o logs/merge_${gene2}.log -N merge_${gene2} \
    6_merge_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF1} ${truncR1}
  echo ${gene2} "denoising is already complete; we will start with merging reads."
elif [[ ! -f "../data/working/${gene2}_7_chimera.RData" ]]; then
  qsub -o logs/chimera_${gene2}.log -N chimera_${gene2} \
    7_chimera_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF1} ${truncR1}
  echo ${gene2} "merging is already complete; we will start with chimera removal."
elif [[ ! -f "../data/working/${gene2}_8_output.RData" ]]; then
  qsub -o logs/output_${gene2}.log -N output_${gene2} \
    8_output_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF1} ${truncR1}
  echo ${gene2} "chimera removal is already complete; we will start with exporting results."
else
  echo "All steps for ${gene2} have already completed"
fi


