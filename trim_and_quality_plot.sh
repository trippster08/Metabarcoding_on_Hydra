# /bin/sh

raw="$1"
gene="$2"
data=${raw}/../


if
  [[ -z "$(ls ${raw}/*.fastq.gz 2>/dev/null | grep fastq.gz)" ]]  
then  
  echo "Correct path to raw read files not entered (*.fastq.gz)"
  exit
fi

if
  [[ -z $2 ]]
then
  echo "Genetic code not entered (should be a number between 1 and 25)"
  exit
fi

if
  [[ ${gene} == "COI" || ${gene} == "COX1" || ${gene} == "CO1" ]]
then
  primerF="/scratch/nmnh_lab/macdonaldk/primers/COImlIntF_spacers.fas" && \
  primerR="/scratch/nmnh_lab/macdonaldk/primers/jgCOIR_spacers.fas"
else
  if
    [[ ${gene} == "12S" || ${gene} == "MiFish" ]]
  then
    primerF="/scratch/nmnh_lab/macdonaldk/primers/MiFish_12SF_spacers.fas" && \
    primerR="/scratch/nmnh_lab/macdonaldk/primers/MiFish_12SR_spacers.fas"
  else
    if
      [[ ${gene} == "18S" ]]
    then
      primerF="/scratch/nmnh_lab/macdonaldk/primers/18SF_spacers.fas" && \
      primerR="/scratch/nmnh_lab/macdonaldk/primers/18SR_spacers.fas"
    else
      echo 'Incorrect reference database. Please enter "COI", "COX1", "CO1", "12S", "MiFish", or "18S"'
      exit
    fi
  fi    
fi

mkdir ${data}/trimmed_reads


qsub -o logs/cutadapt.log \
  -N cutadapt \
trim.job ${data} ${primerF} ${primerR}