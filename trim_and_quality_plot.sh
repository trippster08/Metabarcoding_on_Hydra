# /bin/bash

echo ${@}
echo ${#}
raw=$(readlink -f data/raw)
echo ${raw}
data=$(readlink -f data)
echo ${data}
primerF=$(readlink -f primers/primerF.fas)
primerR=$(readlink -f primers/primerR.fas)
primerFrc=$(readlink -f primers/primerF_RC.fas)
primerRrc=$(readlink -f primers/primerR_RC.fas)
# List of possible primers
possible_primers=(COI 18S MiFish 16S 28SAnth 16Sbac)
echo ${possible_primers[@]}
RC_primers=(MiFish)
echo ${RC_primers[@]}

if ! ls "${raw}"/*.fastq.gz 2>/dev/null | grep -q fastq.gz; then
  echo "No sequences (*.fastq.gz) were found in the raw data directory: ${raw}"
  exit 1
fi

: > primers/primerF.fas
: > primers/primerR.fas
: > primers/primerF_RC.fas
: > primers/primerR_RC.fas

if [ "$#" -eq 1 ]; then
  gene="$1"
  path_to_trimmed=${data}"/working/trimmed_sequences/"
  if [[ ! " ${possible_primers[@]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]]; then
    echo "Error: ${gene} is not a valid primer name. Valid primer names are: ${possible_primers[@]}"
    exit 1
  else
    echo "Illumina run only contains ${gene}. No gene-specific demultiplexing will \
    be performed by cutadapt. If your run contained multiple genes and you expected \
    gene-specific demultiplexing, please kill run using 'qdel -j JOBNUMBER' and include \
    the names of all genes after 'sh trim_and_quality_plot.sh'. Valid gene names are \
    ${possible_primers[@]}."
    if find "${path_to_trimmed}" -name "*.fastq.gz" -print -quit; then
      qsub -o logs/quality.log -N quality \
      2_quality.job
      echo "Trimming is already completed, we moving to the next step: 2_quality.job"
    else
        if [[ " ${RC_primers[*]} " == *" ${gene} "* ]];then
          cp "primers/${gene}F.fas" primers/primerF.fas
          cp "primers/${gene}R.fas" primers/primerR.fas
          cp "primers/${gene}F_RC.fas" primers/primerF_RC.fas
          cp "primers/${gene}R_RC.fas" primers/primerR_RC.fas
          qsub -o logs/cutadapt.log -N cutadapt \
          jobs/1_trim.job ${#} ${@} ${data} ${primerF} ${primerR} ${primerFrc} ${primerRrc}
        else
          cp "primers/${gene}F.fas" primers/primerF.fas
          cp "primers/${gene}R.fas" primers/primerR.fas
          qsub -o logs/cutadapt.log -N cutadapt \
          jobs/1_trim.job ${#} ${@} ${data} ${primerF} ${primerR}
        fi
    fi
  fi
elif [ "$#" -ge 2 ]; then
  path_to_trimmed=${data}"/working/trimmed_sequences/"${1}
  if find "${path_to_trimmed}" -name "*.fastq.gz" -print -quit | grep -q .; then
    echo "Trimming is already completed, we are moving to the next step: 2_quality_multigene.job"
    qsub -o logs/quality_multigene.log -N quality \
    jobs/2_quality_multigene.job ${[@]}
  else
    for gene in "$@"; do
      if [[ ! " ${possible_primers[@]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]]; then
        echo "Error: ${gene} is not a valid primer name. Valid gene names are: \
        ${possible_primers[@]}."
        exit 1
      fi
    done
    if [[ " ${RC_primers[@]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]];then
      cat "primers/${gene}F.fas" >> primers/primerF.fas
      cat "primers/${gene}R.fas" >> primers/primerR.fas
      cat "primers/${gene}F_RC.fas" >> primers/primerF_RC.fas
      cat "primers/${gene}R_RC.fas" >> primers/primerR_RC.fas
      qsub -o logs/cutadapt.log -N cutadapt \
      jobs/1_trim.job ${#} ${@} ${data} ${primerF} ${primerR} ${primerFrc} ${primerRrc}
    else
      cat "primers/${gene}F.fas" >> primers/primerF.fas
      cat "primers/${gene}R.fas" >> primers/primerR.fas
      qsub -o logs/cutadapt.log -N cutadapt \
      jobs/1_trim.job ${#} ${@} ${data} ${primerF} ${primerR}
    fi
  fi  
else
  echo "You must provide at least one gene name following the shell script. Valid \
  gene names are: ${possible_primers[@]}."
  exit 1
fi