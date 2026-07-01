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
  
  -e  <entry='energy_MACE'> entry of the npz file from where the Nmin are
      selected.
  -E  <ending='_epr.out'> ending of the of the files that are excluded in the
      selection of random radicals to compute. In the default case, it excludes
      the radicals that already have an optimization, for example, it excludes
      3 if <file>_003_epr.out already exists.
  -f  <file> file with list of npz files to analyze.
  -j  <job_options='--nice'> joboptions. avoid to add '-n' in this string; that
      variable enters as  the flag -p in of this script.
  -n  <n=1> number of radicals to be selected from each npz file.
  -N  <N=3> number of radicals in the subset that minimizes the entry
      parameter.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
directory='./'
n=1
N=3
entry='energy_MACE'
verbose='false'
priority=''
ending='_epr.out'
restart=''
processors=8
job_options='--nice'
while getopts 'e:E:f:n:N:o:p:j:rvh' flag;
do
  case "${flag}" in
    e) entry=${OPTARG} ;;
    E) ending=${OPTARG} ;;
    f) npz_files=${OPTARG} ;;
    j) job_options=${OPTARG} ;;
    n) n=${OPTARG} ;;
    N) N=${OPTARG} ;;
    o) other_flags=${OPTARG} ;;
    p) processors=${OPTARG} ;;
    r) restart='-R' ;;

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
  fail "There are not npz files defined in the file: $npz_files"
fi

for file in ${all_files[@]}
do
  #Continue computing gvals.
  verbose $file
  existing_n=$(ls "${file%.npz}"_*"$ending" | wc -l 2>/dev/null)
  exclude_indices='[]'
  if [[ $existing_n -gt 0 ]]
  then
    mapfile -t list_of_prev < <(ls "${file%.npz}"_*"$ending" | \
                                   sed "s/$ending//" | \
                                   cut -d _ -f 2 | sed 's/^0*//')
    exclude_indices=$(printf "%d," ${list_of_prev[@]} | sed 's/,$/]/' | \
                      sed 's/^/[/')
  fi
  new_n=$(( n - existing_n ))
  if [[ $new_n -le 0 ]]
  then
    echo "  - $file: already has $existing_n radicals, skipping."
    continue
  fi
  selection=$(myutils nrand_rad $file --n $new_n --Nmin $N \
              --exclude "$exclude_indices") || \
    { warning "Error selecting the radicals for the file: $file"; continue; }
  myutils ext_xyz_from_npz $file --selection "$selection" > /dev/null || failed "extracting the xyz files for the file: $file"

  for i in $(echo "$selection" | tr -d '[],')
  do
    index=$(printf "%03d\n" "$i")
    xyz_file="${file%.npz}_$index.xyz"
    echo "  - $xyz_file"
    info=$(sed -n "2p" $xyz_file | sed "s/;/\n/g")
    charge=$(echo "$info" | grep 'total_charge' | awk '{print $2}')
    multi=$(echo "$info" | grep 'multiplicity' | awk '{print $2}')
    radical=$(echo "$info" | grep 'heavy_atom_missing_idxs' | awk '{print $2}')
    charged_a=$(echo "$info" | grep 'atom_chargerelevant_idx' | awk '{print $2}')
    name=${xyz_file##*/}
    sbatch $job_options -J ${name%.xyz} \
      -n $processors \
      $(myutils gval_workflow -path) -c $charge \
                                     -m $multi \
                                     -f "g$radical,${charged_a}d3" \
                                     -n ${xyz_file%.xyz} \
                                     -b -v -s $restart \
                                     -p $processors $other_flags
  done
done

finish "submitted jobs for the radicals in the npz files in this directory:
  $(pwd)"
