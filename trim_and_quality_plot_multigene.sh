# /bin/bash

raw=$(readlink -f ../data/raw)
cat /dev/null > ../primers/primerF.fas
cat /dev/null > ../primers/primerR.fas
cat /dev/null > ../primers/primerF_RC.fas
cat /dev/null > ../primers/primerR_RC.fas

data=$(readlink -f ../data)
COI=(COI coi CO1 co1 cox1 COX1)
MiFish=(12S MiFish mifish Mifish 12S_mifish 12S_MiFish 12s)
V4=(18S l8s V4 v4)
V4_16S=(16S 16s 16Sbac)
28S=(28S 28s Anth_28S 28S_Anth)

if
  [[ -z "$(ls ${raw}/*.fastq.gz 2>/dev/null | grep fastq.gz)" ]]; then  
  echo "No sequences (*.fastq.gz) were found in the raw data directory: ${raw}"
  exit
fi

# Check if either variable is empty
if [ -z $1 ]; then
  echo "One or both primer sets missing. Please include two of the following three primer sets  "COI", "12S", or "18S"."
  exit 1
fi

# Function to check if a value is in an array

contains() {
  local seeking=$1
  shift
  local arrays=("$@")
  local in=1

  for array in "${arrays[@]}"; do
    for x in "${!array}"; do
      if [[ $x == $seeking ]]; then
        in=0
        break 2
      fi
    done
  done

  return $in
}


found_mifish=0
# Loop through each gene provided as an argument
for gene in "$@"
do
  if contains "$gene" MiFish[@]; then
    found_mifish=1
  fi
done

for gene in "$@"
do
  if [ $found_mifish -eq 1 ]; then
    if contains "$gene" COI[@] MiFish[@] V4[@] V4_16S[@] 28S[@]; then
      cat "../primers/${gene}_primerF.fas" >> ../primers/primerF.fas
      cat "../primers/${gene}_primerR.fas" >> ../primers/primerR.fas
      cat "../primers/${gene}_primerF_RC.fas >> ../primers/primerF_RC.fas
      cat "../primers/${gene}_primerR_RC.fas >> ../primers/primerR_RC.fas
    else
      echo "Primer file for $gene not found."
  fi
done


mkdir -p \
${data}/working/trimmed_sequences/${gene1} \
${data}/working/filtered_sequences/${gene1} \
${data}/working/trimmed_sequences/${gene2} \
${data}/working/filtered_sequences/${gene2} \
${data}/results/${gene1} \
${data}/results/${gene2}

path_to_trimmed=${data}"/working/trimmed_sequences/"${gene1}

if
  [[ -n $(find "${path_to_trimmed}" -name "*.fastq.gz" -print -quit) ]];
then 
  qsub -o logs/quality.log -N quality \
  2_quality_multigene.job ${gene1} ${gene2}
  echo "Trimming is already completed, we moving to the next step: 2_quality_multigene.job"
else
  qsub -o logs/cutadapt.log -N cutadapt \
  1_trim_multigene.job \
  ${data} ${gene1} ${gene2} ${primerF} ${primerR} ${primerFrc} ${primerRrc}
fi