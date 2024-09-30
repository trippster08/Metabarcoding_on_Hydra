# /bin/sh

raw="$1"
gene1="$2"
gene2="$3"
data=${raw}/../


if
  [[ -z "$(ls ${raw}/*.fastq.gz 2>/dev/null | grep fastq.gz)" ]]  
then  
  echo "Correct path to raw read files not entered (*.fastq.gz)"
  exit
fi

# Check if either variable is empty
if [ -z ${gene1} ] || [ -z ${gene2} ]; then
  echo "One or both primers missing. Please include two primer pairs, possibilities include "COI", "COX1", "CO1", "12S", "MiFish", or "18S"."
  exit 1
fi


if
  [[ ${gene1} == "COI" || ${gene1} == "COX1" || ${gene1} == "CO1" ]]
then
  primer1F=${data}"../primers/COImlIntF_spacers.fas" && \
  primer1R=${data}"../primers/jgCOIR_spacers.fas"
else
  if
    [[ ${gene1} == "12S" || ${gene1} == "MiFish" ]]
  then
    primer1F=${data}"../primers/MiFish_12SF_spacers.fas" && \
    primer1R=${data}"../primers/MiFish_12SR_spacers.fas"
  else
    if
      [[ ${gene1} == "18S" ]]
    then
      primer1F=${data}"..//primers/18SF_spacers.fas" && \
      primer1R=${data}"..//primers/18SR_spacers.fas"
    else
      echo 'Incorrect primers referenced. Please enter "COI", "COX1", "CO1", "12S", "MiFish", or "18S"'
      exit
    fi
  fi    
fi

if
  [[ ${gene2} == "COI" || ${gene2} == "COX1" || ${gene2} == "CO1" ]]
then
  primer2F=${data}"../primers/COImlIntF_spacers.fas" && \
  primer2R=${data}"../primers/jgCOIR_spacers.fas"
else
  if
    [[ ${gene2} == "12S" || ${gene2} == "MiFish" ]]
  then
    primer2F=${data}"../primers/MiFish_12SF_spacers.fas" && \
    primer2R=${data}"../primers/MiFish_12SR_spacers.fas"
  else
    if
      [[ ${gene2} == "18S" ]]
    then
      primer2F=${data}"..//primers/18SF_spacers.fas" && \
      primer2R=${data}"..//primers/18SR_spacers.fas"
    else
      echo 'Incorrect primers referenced. Please enter "COI", "COX1", "CO1", "12S", "MiFish", or "18S"'
      exit
    fi
  fi    
fi

mkdir -p ${data}/working/trimmed_reads


qsub -o logs/cutadapt.log \
  -N cutadapt \
trim.job ${raw} ${data} ${primer1F} ${primer1R} ${primer2F} ${primer2R}