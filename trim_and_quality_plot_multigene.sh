# /bin/bash

raw="$1"
gene1="$2"
gene2="$3"
data=${raw}/../
COI=(COI coi CO1 co1 cox1 COX1)
MiFish=(12S MiFish mifish Mifish 12S_mifish 12S_MiFish 12s)
V4=(18S l8s V4 v4)

if
  [[ -z "$(ls ${raw}/*.fastq.gz 2>/dev/null | grep fastq.gz)" ]]; then  
  echo "Correct path to raw read files not entered (*.fastq.gz)"
  exit
fi

# Check if either variable is empty
if [ -z ${gene1} ] || [ -z ${gene2} ]; then
  echo "One or both primer sets missing. Please include two of the following three primer sets  "COI", "12S", or "18S"."
  exit 1
fi

# Function to check if a value is in an array
contains() {
  local array="$1[@]"
  local seeking=$2
  local in=1
  for x in "${!array}"; do
    if [[ $x == $seeking ]]; then
      in=0
      break
    fi
  done
  return $in
}

if (contains COI "$gene1" && contains MiFish "$gene2") || (contains MiFish "$gene1" && contains COI "$gene2"); then
  primerFpath=${data}"../primers/COImlIntF_12SFMiFish_spacers.fas"
  primerRpath=${data}"../primers/jgCOIR_12SRMiFish_spacers.fas"
  primerFrcpath=${data}"../primers/MiFish_12SF_RC.fas"
  primerRrcpath=${data}"../primers/MiFish_12SR_RC.fas"
  primerF=`basename ${primerFpath}`
  primerR=`basename ${primerRpath}`
  primerFrc=`basename ${primerFrcpath}`
  primerRrc=`basename ${primerRrcpath}`
  echo "Primer files set:"
  echo "Forward primer file: `basename ${primerF}`"
  echo "Reverse primer file: $primerR"
  echo "Forward primer RC file: $primerFrc"
  echo "Reverse primer RC file: $primerRrc"
elif (contains V4 "$gene1" && contains MiFish "$gene2") || (contains V4 "$gene2" && contains MiFish "$gene1"); then
  primerFpath=${data}"../primers/18SF_12SFMiFish_spacers.fas"
  primerRpath=${data}"../primers/18SR_12SRMiFish_spacers.fas"
  primerFrcpath=${data}"../primers/MiFish_12SF_RC.fas"
  primerRrcpath=${data}"../primers/MiFish_12SR_RC.fas"
  primerF=`basename ${primerFpath}`
  primerR=`basename ${primerRpath}`
  primerFrc=`basename ${primerFrcpath}`
  primerRrc=`basename ${primerRrcpath}`
  echo "Primer files set:"
  echo "Forward primer file: `basename ${primerF}`"
  echo "Reverse primer file: $primerR"
  echo "Forward primer RC file: $primerFrc"
  echo "Reverse primer RC file: $primerRrc"
elif (contains COI "$gene1" && contains V4 "$gene2") || (contains COI "$gene2" && contains V4 "$gene1"); then
  primerFpath=${data}"../primers/COImlIntF_18SF_spacers.fas"
  primerRpath=${data}"../primers/jgCOIR_18SR_spacers.fas"
  primerF=`basename ${primerFpath}`
  primerR=`basename ${primerRpath}`
  echo "Primer files set:"
  echo "Forward primer file: $primerF"
  echo "Reverse primer file: $primerR"
else
      echo 'No Primer files set. Incorrect primers referenced. Please enter "COI", "12S", or "18S"'
      exit 1   
fi

mkdir -p \
${data}/working/trimmed_sequences/${gene1} \
${data}/working/filtered_sequences/${gene1} \
${data}/working/trimmed_sequences/${gene2} \
${data}/working/filtered_sequences/${gene2} \
${data}/results/${gene1} \
${data}/results/${gene2}

path_to_trimmed=${data}"/working/trimmed_sequences/"${gene1}

if
  [[ -n $(find "$path_to_trimmed" -name "*.fastq.gz" -print -quit) ]];
then 
  qsub -o logs/quality.log -N quality \
  2_quality_multigene.job ${gene1} ${gene2}
  echo "Trimming is already completed, we moving to the next step: 2_quality_multigene.job"
else
  qsub -o logs/cutadapt.log -N cutadapt \
  1_trim_multigene.job \
  ${data} ${gene1} ${gene2} ${primerF} ${primerR} ${primerFrc} ${primerRrc}
fi