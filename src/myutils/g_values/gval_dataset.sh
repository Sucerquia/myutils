#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
This script select n molecules, out of N molecules that minimizes one of the
entries in the npz dile. Then, it goes through the selected molecules,
optimizes them and then computes the g-values and A-matrices (HFC) of the atoms
in the second degree neighborhood of the relevant atoms defined in the npz
files, which are marked as 'heavy_atom_missing_idxs' and
'atom_chargerelevant_idx'.
  
  -d  <directory> location of the dataset containing all the npz files.
  -n  <n=1> number of radicals to be selected from each npz file.
  -N  <N=3> number of radicals in the subset that minimizes the entry
      parameter.
  -e  <entry='energy_MACE'> entry of the npz file from where the Nmin are
      selected.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
directory='./'
n=None
N=None
entry='energy_MACE'
verbose='false'
priority=''
while getopts 'e:f:n:N:pvh' flag;
do
  case "${flag}" in
    f) npz_files=${OPTARG} ;;
    e) entry=${OPTARG} ;;
    n) n=${OPTARG} ;;
    N) N=${OPTARG} ;;
    p) priority='--nice' ;;
  
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
mapfile -t all_files < <(cat $npz_files)


if [[ ${#all_files[@]} -eq 0 ]]
then
  fail "There are not npz files in this directory: $(pwd)"
fi

for file in ${all_files[@]}
do
  verbose $file
  selection=$(myutils nrand_rad $file --n $n --Nmin $N)
  myutils ext_xyz_from_npz $file --selection "$selection" > /dev/null 
  
  for xyz_file in ${file%.npz}_*.xyz
  do
    echo "  - $xyz_file"
    info=$(sed -n "2p" $xyz_file | sed "s/;/\n/g")
    charge=$(echo "$info" | grep 'total_charge' | awk '{print $2}')
    multi=$(echo "$info" | grep 'multiplicity' | awk '{print $2}')
    radical=$(echo "$info" | grep 'heavy_atom_missing_idxs' | awk '{print $2}')
    charged_a=$(echo "$info" | grep 'atom_chargerelevant_idx' | awk '{print $2}')
    name=${xyz_file##*/}
    sbatch --partition=cpu $priority -J ${name%.xyz} \
      $(myutils gval_workflow -path) -c $charge \
                                     -m $multi \
                                     -f "g$radical,${charged_a}d3" \
                                     -n ${xyz_file%.xyz} \
                                     -b -v -s
  done
done

finish "submitted jobs for the radicals in the npz files in this directory:
  $(pwd)"
