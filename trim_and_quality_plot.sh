# /bin/bash
# This prints onto the screen the inputs after "sh trim_and_quality_plot.sh"
#echo ${@}
#echo ${#}
# Set path to raw sequences and data directory
path_to_data=$(cd ./data && pwd)
path_to_raw=${path_to_data}/raw
#echo ${path_to_raw}

#echo ${path_to_data}

# Make a variable to hold all the available primer names
available_primers=()
# Loop through all the files in the primer folder and get the names from
# ecah forward primer file. Then add them to the variable of available primer
# names
for file in primers/*F.fas; do
  primer_file=${file##*/}
  primer_name=${primer_file%F.fas}
  available_primers+=(${primer_name})
done
#echo ${available_primers[@]}

# Do the same for primers with potential readthrough
RC_primers=()
for file in primers/*F_RC.fas; do
  primer_file=${file##*/}
  primer_name=${primer_file%F_RC.fas}
  RC_primers+=(${primer_name})
done
#echo ${RC_primers[@]}

# Look to see if there are fastq.gz files in data/raw
if ! ls "${path_to_raw}"/*.fastq.gz 2>/dev/null | grep -q fastq.gz; then
  echo "No sequences (*.fastq.gz) were found in the raw data directory: ${raw}"
  exit 1
fi

# Create path for trimmed sequences, results, and primers for cutadapt
path_to_trimmed="${path_to_data}/working/trimmed_sequences/"
path_to_results="${path_to_data}/results/"
path_to_primers="primers/active"
# Create directory for active primers
mkdir ${path_to_primers}
# echo ${path_to_trimmed}
# echo ${path_to_results}

# Set path to the 4 potential primer files
PrimerF=${path_to_primers}/primerF.fas
PrimerR=${path_to_primers}/primerR.fas
PrimerFrc=${path_to_primers}/primerF_RC.fas
PrimerRrc=${path_to_primers}/primerR_RC.fas

# Create empty files for the primer sequences to be used by cutadapt. These will
# be filled with sequences below
: > ${PrimerF}
: > ${PrimerR}
: > ${PrimerFrc}
: > ${PrimerRrc}

# This is the loop to submit the job that will perform primer trimming
# (and demultiplexing, if necessary) and create quality plots
# First, see if any variables were submitted?
if [ "$#" -ge 1 ]; then # If variables were submitted
  # Check if all_trimmed = true
  if  [ -d ${path_to_trimmed} ]; then # if the path_to_trimmed does have a folder (but no sequences in them)
    # Set variable for whether all trimming has been completed to true
    all_trimmed=true
    # Check to see if trimming has been completed for all jobs
    # Loop through all genes, determining whether to change all_trimmed to false
    for gene in ${@}; do
      if ! find "${path_to_trimmed}/${gene}" -maxdepth 1 -name "*.fastq.gz" | grep -q .; then
        all_trimmed=false
        break
      fi
    done
    # Check to see if all_trimmed = true (i.e. trimming for all genes is completed)
    if [ ${all_trimmed} = true ]; then # if trimmed sequences for all genes do exist,
    # trimming is complete and the next job, 2_quality.job, is submitted instead 
      qsub -o logs/quality.log -N quality \
      jobs/2_quality.job  ${#} ${@}
      echo "Trimming is already completed, we are moving to the next step: 2_quality.job"
    # If the path_to_trimmed does have a folder (but no sequences in them) Check
    # to see if a cutadapt log and json file exist
    elif find logs/cutadapt* -maxdepth 1 -name '*.json' | grep -q .; then # If
      # both exist, remove both and the trimmed and results folders
      rm -r ${path_to_trimmed} ${path_to_results} logs/cutadapt.log logs/*.json
    elif [  -f logs/cutadapt.log ]; then # If only log exists, remove log and
      # trimmed and results folders
      rm -r ${path_to_trimmed} ${path_to_results} logs/cutadapt.log
    else # If log does not exist, only remove trimmed and results folders
        rm -r ${path_to_trimmed} ${path_to_results}
    fi
  else # If trimmed folder does not exist, make gene-specific trimmed and results folders
    for gene in ${@}; do
      mkdir -p "${path_to_trimmed}${gene}" "${path_to_results}${gene}/plots" \
      "${path_to_results}${gene}/additional_results"
    done
  # Set RC_found to false to start, and only change to true if one of the RC
  # primers is given
    RC_found=false # initialize RC_found outside loop
    # Loop through all the variables (genes) given 
    for gene in "$@"; do
      if [[ ! " ${available_primers[*]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]]; then
        # If one of the variables is not a valid primer, print error and list of primers
        echo "Error: ${gene} is not a valid primer name. Valid gene names are: ${available_primers[@]}."
        exit 1
      fi
      # Check to see if one of the submitted variables is a primer with read-through
      if [[ " ${RC_primers[*]} " == *" ${gene} "* ]];then # If it is, then we need
      # to pass four primers to the cutadapt
        # This adds the sequences from the gene-specific files to the files that will
        # be used by cutadapt. Also set RC_found variable to true
        cat "primers/${gene}F.fas" >> ${PrimerF}
        cat "primers/${gene}R.fas" >> ${PrimerR}
        cat "primers/${gene}F_RC.fas" >> ${PrimerFrc}
        cat "primers/${gene}R_RC.fas" >> ${PrimerRrc}
        RC_found=true
      else # If no RC primers are used, we only need 2 primer files
        # This adds the sequences from the gene-specific files to the files that will
        # be used by cutadapt.
        cat "primers/${gene}F.fas" >> ${PrimerF}
        cat "primers/${gene}R.fas" >> ${PrimerR}
      fi
    done
    #echo ${RC_found}
    # Check to see if RC_found is true
    if [ "$RC_found" = true ]; then # If we used a read-through primer
      # submit job to hydra with primers and rc primers. Also pass number of genes,
      # list of genes, and path to data
      qsub -o logs/cutadapt.log -N cutadapt \
      jobs/1_trim.job ${#} ${@} ${path_to_data} ${PrimerF} ${PrimerR} ${PrimerFrc} ${PrimerRrc}
    else # If no read-through primer
      # submit job to hydra with primers, number of genes, list of genes, and path to data
      qsub -o logs/cutadapt.log -N cutadapt \
      jobs/1_trim.job ${#} ${@} ${path_to_data} ${PrimerF} ${PrimerR}
    fi
  fi
else # If no variables were submitted after shell script
  # Print onto screen that gene names are needed, and gives list of valid primers
  echo "You must provide at least one gene name following the shell script. Valid \
  gene names are: ${available_primers[@]}."
  exit 1
fi