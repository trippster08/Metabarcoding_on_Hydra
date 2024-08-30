# /bin/sh

raw=$"1"

# Please change "Path_to_forward_primer_including_primer_file" to the path of
# your forward primer, and "Path_to_reverse_primer_including_primer_file" to 
# the path of your reverse primer. Don't use relative paths here, for simplicity
primerF="Path_to_forward_primer_including_primer_file"
primerR="Path_to_reverse_primer_including_primer_file"
# Please change 



qsub -o logs/${raw}_cutadapt.log \
  -N ${raw}_fastp \
cutadapt.job ${primerF} ${primerR}