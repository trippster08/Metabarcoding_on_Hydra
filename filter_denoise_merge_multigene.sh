# /bin/bash
trimmed1="$1"
trimmed2="$2"
truncF1="$3"
truncR1="$4"
truncF2="$5"
truncR2="$6"

echo ${trimmed1}
echo ${trimmed2}


if
  [[ -z "$(ls ${trimmed1}/*.fastq.gz 2>/dev/null | grep fastq)" ]] || \
  [[ -z "$(ls ${trimmed2}/*.fastq.gz 2>/dev/null | grep fastq)" ]]  
then  
  echo "Correct path to at least one of the directories holding your trimmed read files not entered"
  exit
fi

if
  [[ ! ${truncF1} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncR1} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncF2} =~ ^[0-9]+$ ]] || \
  [[ ! ${truncR2} =~ ^[0-9]+$ ]]
then
  echo "Length to truncation values for one of your genes is not entered correctly. Please enter '<path_to_trimmed_gene1_reads> <path_to_trimmed_gene2_reads> <R1_truncation_value_for_gene1> <R2_truncation_value_for_gene1>  <R1_truncation_value_for_gene2> <R2_truncation_value_for_gene2>'."
  exit
fi

gene1=`basename ${trimmed1}`
gene2=`basename ${trimmed2}`

mkdir -p \
${trimmed1}/../../../results/${gene1} \
${trimmed1}/../../../results/${gene2}

#echo ${gene1}
#echo ${gene2}
#echo ${trimmed1}
#echo ${trimmed2}
#echo ${truncF1}
#echo ${truncR1}
#echo ${truncF2}
#echo ${truncR2}

qsub -o logs/denoise.log \
-N denoise \
filter_denoise_merge_multigene.job ${gene1} ${gene2} ${trimmed1} ${trimmed2} \
${truncF1} ${truncR1} ${truncF2} ${truncR2}