# /bin/sh

raw="$1"
gene="$2"
data=${raw}/../
COI=(COI coi CO1 co1 cox1 COX1)
MiFish=(12S MiFish mifish Mifish 12S_mifish 12S_MiFish 12s)
V4=(18S l8s V4 v4)

if
  [[ -z "$(ls ${raw}/*.fastq.gz 2>/dev/null | grep fastq.gz)" ]]; then  
  echo "Correct path to raw read files not entered (*.fastq.gz)"
  exit
fi

if
  [[ -z ${gene} ]]
then
  echo "Primer set not given (Please enter "COI", "COX1", "CO1", "12S", "MiFish", or "18S" after path to raw reads)"
  exit
fi

if [[ " ${COI[@]} " =~ " ${gene} " ]]; then
  primerFpath=${data}"../primers/COImlIntF_spacers.fas"
  primerRpath=${data}"../primers/jgCOIR_spacers.fas"
  primerF=`basename ${primerFpath}`
  primerR=`basename ${primerRpath}`
elif [[ " ${MiFish[@]} " =~ " ${gene} " ]]; then
  primerFpath=${data}"../primers/MiFish_12SF_spacers.fas"
  primerRpath=${data}"../primers/MiFish_12SR_spacers.fas"
  primerFrcpath=${data}"../primers/MiFish_12SF_RC.fas"
  primerRrcpath=${data}"../primers/MiFish_12SR_RC.fas"
  primerF=`basename ${primerFpath}`
  primerR=`basename ${primerRpath}`
  primerFrc=`basename ${primerFrcpath}`
  primerRrc=`basename ${primerRrcpath}`
elif [[ " ${V4[@]} " =~ " ${gene} " ]]; then
  primerFpath=${data}"../primers/18SF_spacers.fas"
  primerRpath=${data}"../primers/18SR_spacers.fas"
  primerF=`basename ${primerFpath}`
  primerR=`basename ${primerRpath}`
else
  echo 'Incorrect primer set given. Please enter "COI", "12S", or "18S"'
  exit 1
fi   

# Create all the subdirectories we will use
mkdir -p ../data/working/trimmed_sequences ../data/results

path_to_trimmed=${data}"/working/trimmed_sequences/"${gene1}

if
  [[ -n $(find "$path_to_trimmed" -name "*.fastq.gz" -print -quit) ]];
then
  qsub -o logs/quality.log -N quality \
  2_quality.job
  echo "Trimming is already completed, we moving to the next step: 2_quality.job"
else
  qsub -o logs/cutadapt.log -N cutadapt \
  1_trim.job ${data} ${primerF} ${primerR} ${primerFrc} ${primerRrc}
fi