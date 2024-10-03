# /bin/sh

raw="$1"
primer="$2"
data=${raw}/../

wget https://github.com/trippster08/Metabarcoding_on_Hydra/archive/refs/heads/main.zip

unzip main.zip
mv Metabarcoding_on_Hydra-main/* .
mv primers ..
rm -r Metabarcoding_on_Hydra-main

if
  [[ -z "$(ls ${raw}/*.fastq.gz 2>/dev/null | grep fastq.gz)" ]]  
then  
  echo "Correct path to raw read files not entered (*.fastq.gz)"
  exit
fi

if
  [[ -z $2 ]]
then
  echo "Primer set not given (Please enter "COI", "COX1", "CO1", "12S", "MiFish", or "18S" after path to raw reads)"
  exit
fi

if
  [[ ${gene} == "COI" || ${gene} == "COX1" || ${gene} == "CO1" ]]
then
  primerF=${data}"../primers/COImlIntF_spacers.fas" && \
  primerR=${data}"../primers/jgCOIR_spacers.fas"
else
  if
    [[ ${gene} == "12S" || ${gene} == "MiFish" ]]
  then
    primerF=${data}"../primers/MiFish_12SF_spacers.fas" && \
    primerR=${data}"../primers/MiFish_12SR_spacers.fas"
  else
    if
      [[ ${gene} == "18S" ]]
    then
      primerF=${data}"../primers/18SF_spacers.fas" && \
      primerR=${data}"../primers/18SR_spacers.fas"
    else
      echo 'Incorrrect primer set given. Please enter "COI", "COX1", "CO1", "12S", "MiFish", or "18S"'
      exit
    fi
  fi    
fi

mkdir -p ${data}/working/trimmed_reads


qsub -o logs/cutadapt.log \
  -N cutadapt \
trim_and_quality_plot.job ${data} ${primerF} ${primerR}