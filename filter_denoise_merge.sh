# /bin/sh
trimmed="$1"
truncF="$2"
truncR="$3"
data=${trimmed}/../../






qsub -o logs/denoise.log \
  -N denoise \
trim.job ${trimmed} ${truncF} ${truncR} ${data}