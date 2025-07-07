# /bin/bash

#echo ${@}
#echo ${#}
raw=$(readlink -f data/raw)
#echo ${raw}
data=$(readlink -f data)
#echo ${data}
primerF=$(readlink -f primers/primerF.fas)
primerR=$(readlink -f primers/primerR.fas)
primerFrc=$(readlink -f primers/primerF_RC.fas)
primerRrc=$(readlink -f primers/primerR_RC.fas)
# List of possible primers
possible_primers=(COI 18S MiFish 16S 28SAnth 16Sbac)
#echo ${possible_primers[@]}
RC_primers=(MiFish)
#echo ${RC_primers[@]}

if ! ls "${raw}"/*.fastq.gz 2>/dev/null | grep -q fastq.gz; then
  echo "No sequences (*.fastq.gz) were found in the raw data directory: ${raw}"
  exit 1
fi

: > primers/primerF.fas
: > primers/primerR.fas
: > primers/primerF_RC.fas
: > primers/primerR_RC.fas


if [ "$#" -ge 1 ]; then
  path_to_trimmed="${data}/working/trimmed_sequences"
  path_to_results="${data}/results"
  if [ -d ${path_to_trimmed} ]; then
    if find "${path_to_trimmed}" -name "*.fastq.gz" | grep -q .; then
      qsub -o logs/quality.log -N quality \
      jobs/2_quality.job
      echo "Trimming is already completed, we are moving to the next step: 2_quality.job"
    else
      if find logs/cutadapt.log -maxdepth 1 -name '*.json' | grep -q .; then 
        rm -r ${path_to_trimmed} logs/cutadapt.log logs/*.json
      elif [  -f logs/cutadapt.log ]; then
        rm -r ${path_to_trimmed} logs/cutadapt.log
      else
          rm -r ${path_to_trimmed}
      fi
    fi
  else
    mkdir -p ${path_to_trimmed} ${path_to_results}
  fi

  RC_found=false # initialize RC_found outside loop

  for gene in "$@"; do
    if [[ ! " ${possible_primers[@]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]]; then
      echo "Error: ${gene} is not a valid primer name. Valid gene names are: ${possible_primers[@]}."
      exit 1
    else
      mkdir -p ${path_to_trimmed}/${gene}
    fi
    if [[ " ${RC_primers[*]} " == *" ${gene} "* ]];then
      cat "primers/${gene}F.fas" >> primers/primerF.fas
      cat "primers/${gene}R.fas" >> primers/primerR.fas
      cat "primers/${gene}F_RC.fas" >> primers/primerF_RC.fas
      cat "primers/${gene}R_RC.fas" >> primers/primerR_RC.fas
      RC_found=true
    else
      cat "primers/${gene}F.fas" >> primers/primerF.fas
      cat "primers/${gene}R.fas" >> primers/primerR.fas
    fi
  done
  #echo ${RC_found}

  if [ "$RC_found" = true ]; then
    qsub -o logs/cutadapt.log -N cutadapt \
    jobs/1_trim.job ${#} ${@} ${data} ${primerF} ${primerR} ${primerFrc} ${primerRrc}
  else
    qsub -o logs/cutadapt.log -N cutadapt \
    jobs/1_trim.job ${#} ${@} ${data} ${primerF} ${primerR}
  fi
else
  echo "You must provide at least one gene name following the shell script. Valid \
  gene names are: ${possible_primers[@]}."
  exit 1
fi