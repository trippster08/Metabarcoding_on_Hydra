# /bin/bash

cat /dev/null > ../primers/primerF_RC.fas
cat /dev/null > ../primers/primerR_RC.fas
echo $1
echo $2

# List of possible primers
possible_primers=(COI 18S MiFish 16S 28SAnth 16Sbac)
echo ${possible_primers[@]}
RC_primers=(MiFish)
echo ${RC_primers[@]}
# Check if the number of arguments is zero or one
if [ "$#" -lt 2 ]; then
  echo "You must provide more than one gene name."
  exit 1
fi

if [ -f ../primers/primerF.fas ];then
  > ../primers/primerF.fas
  > ../primers/primerR.fas
else
  cat /dev/null > ../primers/primerF.fas
  cat /dev/null > ../primers/primerR.fas
fi

if [ -f ../primers/primerF_RC.fas ];then
  > ../primers/primerF_RC.fas
  > ../primers/primerR._RCfas
fi

create_RC_files=false

# Loop through each gene name provided as argument
for gene in "$@"; do
  if [[ ! " ${possible_primers[@]} " =~ " ${gene} " ]]; then
    echo "Error: ${gene} is not a valid primer name. Valid primer names are: ${possible_primers[@]}"
    exit 1
  else
    if [[ " ${RC_primers[@]} " =~ " ${gene} " ]];then
      cat "../primers/${gene}F.fas" >> ../primers/primerF.fas
      cat "../primers/${gene}R.fas" >> ../primers/primerR.fas
      cat "../primers/${gene}F_RC.fas" >> ../primers/primerF_RC.fas
      cat "../primers/${gene}R_RC.fas" >> ../primers/primerR_RC.fas
    else
      cat "../primers/${gene}F.fas" >> ../primers/primerF.fas
      cat "../primers/${gene}R.fas" >> ../primers/primerR.fas
    fi
  fi
done



