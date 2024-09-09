# /bin/sh
trimmed="$1"
truncF="$2"
truncR="$3"

if
  [[ -z "$(ls ${trimmed}/*.fastq.gz 2>/dev/null | grep fastq.gz)" ]]  
then  
  echo "Correct path to trimmed read files not entered (*.fastq.gz)"
  exit
fi

if
  [[ -z $2 ]]
then
  echo "Length to truncate forward reads not entered (should be a number > 1)"
  exit
fi
if
  [[ -z $3 ]]
then
  echo "Length to truncate reverse reads not entered (should be a number > 1)"
  exit
fi


qsub -o logs/denoise.log \
  -N denoise \
filter_denoise_merge.job ${trimmed} ${truncF} ${truncR}