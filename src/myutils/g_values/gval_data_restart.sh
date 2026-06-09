#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
takes a list of xyz files and restarts the gval workflow for those files. This
is useful when the previous run failed or was interrupted.
  
  -f  <file> file with list of xyz files to restart in the current directory.
  -p  <n processors> number of processors.
  -J  <job options> additional job options, e.g., \"--partition=cpu\".
  -P  use for preemptionable jobs.
  -o  add '-R' to restart, '-o' to avoid optimization or '-b' to avoid
      frequencies.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
job_options=""
preemption=""
while getopts 'p:f:J:o:Pvh' flag;
do
  case "${flag}" in
    f) xyz_files=${OPTARG} ;;
    p) processors=${OPTARG} ;;
    J) job_options=${OPTARG} ;;
    o) others_flags=${OPTARG} ;;
    P) preemption='-P' ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" submit_dataset $verbose

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# ---- BODY -------------------------------------------------------------------

# finds all the npz files in the directory and subdirectories
mapfile -t all_files < <(cat $xyz_files)

for xyz_file in ${all_files[@]}
do
  echo "  - $xyz_file"
  info=$(sed -n "2p" $xyz_file | sed "s/;/\n/g")
  charge=$(echo "$info" | grep 'total_charge' | awk '{print $2}')
  multi=$(echo "$info" | grep 'multiplicity' | awk '{print $2}')
  radical=$(echo "$info" | grep 'heavy_atom_missing_idxs' | awk '{print $2}')
  charged_a=$(echo "$info" | grep 'atom_chargerelevant_idx' | awk '{print $2}')
  name=${xyz_file##*/}

  if [[ "$preemption" == "-P" ]]
  then
    echo $job_options | grep -q qos || job_options="$job_options --qos=low"
  fi

  sbatch -n $processors $priority -J ${name%.xyz}  $job_options \
    $(myutils gval_workflow -path) -c $charge \
                                   -m $multi \
                                   -f "g$radical,${charged_a}d3" \
                                   -n ${xyz_file%.xyz} \
                                   -b -v -s -p $processors $preemption $other_flags
done

finish "submitted jobs for the radicals in the npz files in this directory:
  $(pwd)"
