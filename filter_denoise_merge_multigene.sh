 # /bin/bash
trimmed="$1"
gene1="$2"
gene2="$3"
truncF1="$4"
truncR1="$5"
truncF2="$6"
truncR2="$7"
COI=(COI coi CO1 co1 cox1 COX1)
MiFish=(12S MiFish mifish Mifish 12S_mifish 12S_MiFish 12s)
V4=(18S l8s V4 v4)

echo "Your trimmed reads are in this directory:" ${trimmed}

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


if
  [[ -z "$(ls ${trimmed}/${gene1}/*.fastq.gz 2>/dev/null | grep fastq)" ]] || \
  [[ -z "$(ls ${trimmed}/${gene2}/*.fastq.gz 2>/dev/null | grep fastq)" ]]  
then  
  echo "Correct path to at least one of the directories holding your trimmed read files not entered"
  exit
fi

if
  [[ ! ${truncF1} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncR1} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncF2} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncR2} =~ ^[0-9]+$ ]]
then
  echo "Length to truncation values for one of your genes is not entered correctly. Please enter <path_to_trimmed_sequences> <gene1> <gene2> <R1_truncation_value_for_gene1> <R2_truncation_value_for_gene1>  <R1_truncation_value_for_gene2> <R2_truncation_value_for_gene2>'."
  exit
fi

mkdir -p \
${trimmed}/../../results/${gene1} \
${trimmed}/../../results/${gene2}

#echo ${gene1}
#echo ${gene2}
#echo ${trimmed}
#echo ${truncF1}
#echo ${truncR1}
#echo ${truncF2}
#echo ${truncR2}

if [[ ! -f "../data/working/${gene1}_3_filter.RData" ]]; then
  qsub -o logs/${gene1}_filter.log -N ${gene1}_filter \
    filter_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
elif [[ ! -f "../data/working/${gene1}_4_error.RData" ]]; then
  qsub -o logs/${gene1}_error.log -N ${gene1}_error \
    error_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
elif [[ ! -f "../data/working/${gene1}_5_denoise.RData" ]]; then
  qsub -o logs/${gene1}_denoise.log -N ${gene1}_denoise \
    denoise_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
elif [[ ! -f "../data/working/${gene1}_6_merge.RData" ]]; then
  qsub -o logs/${gene1}_merge.log -N ${gene1}_merge \
    merge_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
elif [[ -f "../data/working/${gene1}_7_chimera.RData" ]]; then
  qsub -o logs/${gene1}_chimera.log -N ${gene1}_chimera \
    chimera_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
elif [[ -f "../data/working/${gene1}_8_output.RData" ]]; then
  qsub -o logs/${gene1}_output.log -N ${gene1}_output \
  output_multigene.job ${gene1} ${trimmed}/${gene1} ${truncF1} ${truncR1}
else
  echo "All steps for ${gene1} have already completed"
fi

if [[ ! -f "../data/working/${gene2}_3_filter.RData" ]]; then
  qsub -o logs/${gene2}_filter.log -N ${gene2}_filter \
    filter_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF2} ${truncR2}
elif [[ ! -f "../data/working/${gene2}_4_error.RData" ]]; then
  qsub -o logs/${gene2}_error.log -N ${gene2}_error \
    error_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF2} ${truncR2}
elif [[ ! -f "../data/working/${gene2}_5_denoise.RData" ]]; then
  qsub -o logs/${gene2}_denoise.log -N ${gene2}_denoise \
    denoise_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF2} ${truncR2}
elif [[ ! -f "../data/working/${gene2}_6_merge.RData" ]]; then
  qsub -o logs/${gene2}_merge.log -N ${gene2}_merge \
    merge_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF2} ${truncR2}
elif [[ -f "../data/working/${gene2}_7_chimera.RData" ]]; then
  qsub -o logs/${gene2}_chimera.log -N ${gene2}_chimera \
    chimera_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF2} ${truncR2}
elif [[ -f "../data/working/${gene2}_8_output.RData" ]]; then
  qsub -o logs/${gene2}_output.log -N ${gene2}_output \
  output_multigene.job ${gene2} ${trimmed}/${gene2} ${truncF2} ${truncR2}
else
  echo "All steps for ${gene2} have already completed"
fi


