# /bin/sh
trimmed="$1"
truncF="$2"
truncR="$3"

if
  [[ -z "$(ls ${trimmed}/*.fastq.gz 2>/dev/null | grep fastq)" ]]  
then  
  echo "Correct path to trimmed read files not entered (*.fastq.gz)"
  exit
fi

if
  [[ ! ${truncF} =~ ^[0-9]+$ ]] || \
    [[ ! ${truncR} =~ ^[0-9]+$ ]]
then
  echo "At least one of the lengths-to-truncate values for (values should be a number > 1). Please enter '<path_to_trimmed_reads> <R1_truncation_value> <R2_truncation_value>" 
  exit
fi

qsub -o logs/denoise.log \
  -N denoise \
filter_denoise_merge.job ${trimmed} ${truncF} ${truncR}