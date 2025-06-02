#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Use this template to create your scripts with a standard structure
  
  -f  <file> npz file containing radicals information.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
npz_file=''
verbose='false'
while getopts 'f:vh' flag;
do
  case "${flag}" in
    f) file=${OPTARG} ;;
  
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

myutils ext_xyz_from_npz $file > /dev/null

for xyz_file in ${file%.npz}_*.xyz
do
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
                                       -v
done

finish "message to finish"
