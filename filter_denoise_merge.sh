# /bin/sh
trimmed="realpath ../data/working/raw"
truncF="$2"
truncR="$3"

if
  [[ -z "$(ls ${trimmed}/*.fastq.gz 2>/dev/null | grep fastq)" ]]  
then  
  echo "No sequences (*.fastq.gz) were found in the trimmed data directory: $(realpath ../data/working/trimmed_sequences)"
  exit
fi

if
  [[ ! ${truncF} =~ ^[0-9]+$ ]] || \
    [[ ! ${truncR} =~ ^[0-9]+$ ]]
then
  echo "At least one of the lengths-to-truncate values for (values should be a number > 1). Please enter '<path_to_trimmed_reads> <R1_truncation_value> <R2_truncation_value>" 
  exit
fi

if [[ ! -f "../data/working/3_filter.RData" ]]; then
  qsub -o logs/filter.log -N filter \
    3_filter.job ${trimmed} ${truncF} ${truncR}
elif [[ ! -f "../data/working/4_error.RData" ]]; then
  qsub -o logs/error.log -N error \
    4_error.job ${trimmed} ${truncF} ${truncR}
  echo "Filtering is already complete; we will start with error modelling."
elif [[ ! -f "../data/working/5_denoise.RData" ]]; then
  qsub -o logs/denoise.log -N denoise \
    5_denoise.job ${trimmed} ${truncF} ${truncR}
  echo "Error modelling is already complete; we will start with denoising."
elif [[ ! -f "../data/working/6_merge.RData" ]]; then
  qsub -o logs/merge.log -N merge \
    6_merge.job ${trimmed} ${truncF} ${truncR}
  echo "Denoising is already complete; we will start with merging reads."
elif [[ ! -f "../data/working/7_chimera.RData" ]]; then
  qsub -o logs/chimera.log -N chimera \
    7_chimera.job ${trimmed} ${truncF} ${truncR}
  echo "Merging is already complete; we will start with chimera removal."
elif [[ ! -f "../data/working/8_output.RData" ]]; then
  qsub -o logs/output.log -N output \
  8_output.job ${trimmed} ${truncF} ${truncR}
  echo "Chimera removal is already complete; we will start with exporting results."
else
  echo "All steps have already completed"
fi