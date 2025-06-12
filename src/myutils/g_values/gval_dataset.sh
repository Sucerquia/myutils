#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
This script goes through a set of molecules stored in npz files, optimizes them
and then computes the g-values and A-matrices (HFC) of the atoms in the second
degree neighborhood of the relevant atoms defined in the npz files, which are
marked as 'heavy_atom_missing_idxs' and 'atom_chargerelevant_idx'.
  
  -d  <directory> location of the dataset containing all the npz files.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
directory='./'
verbose='false'
while getopts 'd:vh' flag;
do
  case "${flag}" in
    d) directory=${OPTARG} ;;
  
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
cd $directory
# finds all the npz files in the directory and subdirectories
mapfile -t all_files < <(find . -name '*.npz')

if [[ ${#all_files[@]} -eq 0 ]]
then
  fail "There are not npz files in this directory: $(pwd)"
fi

for file in ${all_files[@]}
do
  verbose $file
  myutils ext_xyz_from_npz $file > /dev/null
  
  for xyz_file in ${file%.npz}_*.xyz
  do
    echo "\ \ - $xyz_file"
    info=$(sed -n "2p" $xyz_file | sed "s/;/\n/g")
    charge=$(echo "$info" | grep 'total_charge' | awk '{print $2}')
    multi=$(echo "$info" | grep 'multiplicity' | awk '{print $2}')
    radical=$(echo "$info" | grep 'heavy_atom_missing_idxs' | awk '{print $2}')
    charged_a=$(echo "$info" | grep 'atom_chargerelevant_idx' | awk '{print $2}')
    
    sbatch --partition=cpu-single -J ${xyz_file%.xyz} \
      $(myutils gval_workflow -path) -c $charge \
                                     -m $multi \
                                     -f "g$radical,${charged_a}d3" \
                                     -n ${xyz_file%.xyz} \
                                     -b -v -s
  done
done

finish "message to finish"
