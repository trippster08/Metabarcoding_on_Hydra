# /bin/sh

for x in ../../dataset*; do 
  sample=`basename ${x}`
  echo ${sample}

qsub -o logs/${sample}_tax_assign.log \
-N ${sample}_tax_assign \
tax_assign.job ${sample} 
  sleep 0.1
done