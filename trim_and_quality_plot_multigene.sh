# /bin/sh

raw="$1"
gene1="$2"
gene2="$3"
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

# Check if either variable is empty
if [ -z ${gene1} ] || [ -z ${gene2} ]; then
  echo "One or both primer sets missing. Please include two of the following three primer sets  "COI", "12S", or "18S"."
  exit 1
fi


if
  [[ ${gene1} == "COI" || ${gene1} == "12S" && ${gene2} == "COI" || ${gene2} == "12S" ]]
then
  primerF=${data}"../primers/COImlIntF_12SFMiFish_spacers.fas" && \
  primerR=${data}"../primers/jgCOIR_12SRMiFish_spacers.fas" && \
  primerFrc=${data}"../primers/MiFish_12SF_RC_spacers.fas" && \
  primerRrc=${data}"../primers/MiFish_12SR_RC_spacers.fas"
else
  if
  [[ ${gene1} == "18S" || ${gene1} == "12S" && ${gene2} == "18S" || ${gene2} == "12S" ]]
  then
  primerF=${data}"../primers/18SF_12SFMiFish_spacers.fas" && \
  primerR=${data}"../primers/18SR_12SRMiFish_spacers.fas" && \
  primerFrc=${data}"../primers/MiFish_12SF_RC_spacers.fas" && \
  primerRrc=${data}"../primers/MiFish_12SR_RC_spacers.fas"
  else
  [[ ${gene1} == "COI" || ${gene1} == "18S" && ${gene2} == "COI" || ${gene2} == "18S" ]]
then
  primerF=${data}"../primers/COImlIntF_18SF_spacers.fas" && \
  primerR=${data}"../primers/jgCOIR_18SR_spacers.fas" &&
    else
      echo 'Incorrect primers referenced. Please enter "COI", "12S", or "18S"'
      exit
    fi
  fi    
fi

mkdir -p \
${data}/working/trimmed_reads/${gene1} \
${data}/working/filtered_reads/${gene1} \
${data}/working/trimmed_reads/${gene2} \
${data}/working/filtered_reads/${gene2} \

qsub -o logs/cutadapt.log \
  -N cutadapt \
trim_and_quality_plot_multigene.job ${raw} ${data} ${primerF} ${primerR} ${primerFrc} ${primerRrc} ${gene1} ${gene2}