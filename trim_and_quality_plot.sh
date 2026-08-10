#!/bin/bash

# Set path to raw sequences and data directory
path_to_data=$(realpath ./data)
path_to_raw=${path_to_data}/raw
#echo ${path_to_raw}

#echo ${path_to_data}

# Make a variable to hold all the available primer names
available_primers=()
# Loop through all the files in the primer folder and get the names from
# ecah forward primer file. Then add them to the variable of available primer
# names
for file in primers/*-F.fas; do
  primer_file=${file##*/}
  primer_name=${primer_file%-F.fas}
  available_primers+=(${primer_name})
done
#echo ${available_primers[@]}

# Look to see if there are fastq.gz files in data/raw. 
if ! compgen -G "${path_to_raw}"/*.fastq.gz > /dev/null; then
  echo "No sequences (*.fastq.gz) were found in the raw data directory: ${path_to_raw}."
  echo "Looking for fastq.gz files in other project directories."

  # If no fastq.gz in data/raw, search for the fastq.gz files in other project directories
  filelist="/tmp/fastq_list_$$"
  found_fastq=0
  # Make the raw directory path relative to current directory (so find can exclude it)
  raw_rel="./$(realpath --relative-to="." "$path_to_raw")"

  find -L . \
    -mindepth 2 \
    -type f \
    -name '*.fastq.gz' \
    ! -iname '*undetermined*' \
    ! -path "${raw_rel}/*" \
    -print0 > "$filelist"
  # If no fastq.gz files are found
  if [ ! -s "$filelist" ]; then
    echo "No fastq.gz files were found in any project directories. Please add your raw sequences to the data/raw directory and rerun the script."
    rm -f "$filelist"
    exit 1
  fi 
  
  while IFS= read -r -d '' file; do
    if [ "$found_fastq" -eq 0 ]; then
      # Get directory of the first found fastq.gz file
      first_dir=$(dirname -- "$file")
      echo "Found fastq.gz files in the following locations: ${first_dir}"
      echo "Copying fastq.gz files to ${path_to_raw}."
      found_fastq=1
    fi
    cp -n -- "$file" "${path_to_raw}/"
  done < "$filelist"
  rm -f "$filelist"
  echo "Finished copying fastq.gz files to ${path_to_raw}."
fi  

# count number of files in data/raw and print to screen
count=$(ls -1 "${path_to_raw}"/*.fastq.gz 2>/dev/null | wc -l)
echo "Total fastq.gz files in ${path_to_raw}: $count"

# Create path for demulitplexed sequences, trimmed sequences,
# results, and primers for cutadapt
path_to_demultiplexed="${path_to_data}/raw/demultiplexed_sequences"
path_to_trimmed="${path_to_data}/working/trimmed_sequences"
path_to_results="${path_to_data}/results"
path_to_primers="primers/active"
# Create directory for active primers, but it only makes it if it does not
# already exist, it won't give an error if it does exist
mkdir -p "${path_to_primers}"
# echo ${path_to_trimmed}
# echo ${path_to_results}
# Remove any primer files that are currently in the activep primer folder before
# adding the new primer files
rm -f "${path_to_primers}"/*.fas
# Set path to the two 5' primer files
PrimerF=${path_to_primers}/primerF.fas
PrimerR=${path_to_primers}/primerR.fas
# Create empty files for the 5' primer sequences to be used by cutadapt. These will
# be filled with sequences below
: > "${PrimerF}"
: > "${PrimerR}"

# This is the loop to submit the job that will perform primer trimming
# (and demultiplexing, if necessary) and create quality plots
# First, see if any variables were submitted?
if [ "$#" -ge 1 ]; then # If variables were submitted
  # Check if all_trimmed = true
  if  [ -d "${path_to_trimmed}" ]; then # if the path_to_trimmed does have a folder
    # Set variable for whether all trimming has been completed to true
    all_trimmed=true
    # Check to see if trimming has been completed for all genes
    # Loop through all genes to see if there are sequences in the trimmed_sequences
    # folder. If not we change all_trimmed to false
    for gene in "$@"; do
      if ! find "${path_to_trimmed}/${gene}" -maxdepth 1 -name "*.fastq.gz" | grep -q .; then
        all_trimmed=false
        break
      fi
    done
    # Check to see if all_trimmed = true (i.e. trimming for all genes is completed)
    if [ "${all_trimmed}" = true ]; then # if trimmed sequences for all genes do exist,
    # trimming is complete and the next job, 2_quality.job, is submitted instead 
      qsub -o logs/quality.log -N quality \
      jobs/2_quality.job  ${#} "$@"
      echo "Trimming is already completed, we are moving to the next step: 2_quality.job"
      exit 0
    # If the path_to_trimmed does have a folder (but at least one gene has no sequences in it),
    # check to see if a cutadapt log exists.
    elif [  -f logs/cutadapt.log ]; then # If the log exists, remove log and
      # trimmed and results folders
      rm -rf "${path_to_trimmed}" "${path_to_results}" logs/cutadapt.log
    else # If log does not exist, only remove trimmed and results folders
        rm -rf "${path_to_trimmed}" "${path_to_results}"
    fi
  fi
else # If no variables were submitted after shell script
  # Print onto screen that gene names are needed, and gives list of valid primers
  echo "You must provide at least one gene name following the shell script. Valid \
  gene names are: ${available_primers[@]}."
  exit 1
fi
for gene in "$@"; do
  mkdir -p "${path_to_trimmed}/${gene}" \
  "${path_to_results}/${gene}/plots" \
  "${path_to_results}/${gene}/additional_results" \
  "${path_to_demultiplexed}/${gene}"
done
# We need to set up the primer files for cutadapt before submitting the cutadapt job.
# Loop through all the variables (genes) given 
for gene in "$@"; do
  if [[ ! " ${available_primers[*]} " =~ (^|[[:space:]])${gene}([[:space:]]|$) ]]; then
    # If one of the variables is not a valid primer, print error and list of primers
    echo "Error: ${gene} is not a valid primer name. Valid gene names are: ${available_primers[@]}."
    exit 1
  fi
  # This adds the sequences from the gene-specific files to the files that will
  # be used by cutadapt. Also set RC_found variable to true
  cat "primers/${gene}-F.fas" >> "${PrimerF}"
  cat "primers/${gene}-R.fas" >> "${PrimerR}"

  # Check to see if each of the submitted variables is a primer with read-through
  if [[ -f "primers/${gene}-F_RC.fas" && \
    -f "primers/${gene}-R_RC.fas" ]]; then 
    # If it is, then we need to copy the RC primers into the active primer folder
    # This copies each gene-specific 3' primer file to the active primers folder, and sets the path to
    # each of the gene-specific 3' primer files to be used by cutadapt
    cp "primers/${gene}-F_RC.fas" "${path_to_primers}/${gene}-F_RC.fas"
    cp "primers/${gene}-R_RC.fas" "${path_to_primers}/${gene}-R_RC.fas"
  fi
done
# Submit the job, which will run cutadapt and create quality plots. We are also passing
# the path to the data, the 5' primer files, and the path to the active primers folder,
# as well as the number of variables (genes) and the variables (gene names) that were
# submitted by the user
  qsub -o logs/cutadapt.log -N cutadapt \
  jobs/1_trim.job "${#}" "$@" "${path_to_data}" "${path_to_primers}" "${PrimerF}" "${PrimerR}"

