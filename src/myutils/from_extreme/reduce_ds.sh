#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool finds the dofs and removes all of them that are repeated. such that
the trajectory is reduced without loosing information.

    -d    <directory path = ./>. path to the data set of the peptides.
    -f    <xyz files pattern>. The code will look for *pattern*.xyz
    -p    <peptide>. Alternative to xyz files, this code can extract the peptides
          from <peptide>-optext.log

    -h    prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- general setup ---------------------------------------------------------
pep=''
xyz_pattern=''

ds_dir="./"
while getopts 'd:f:p:s:h' flag;
do
  case "${flag}" in
    d) ds_dir=${OPTARG} ;;
    f) xyz_pattern=${OPTARG} ;;
    p) pep=${OPTARG} ;;
    s) subdir=${OPTARG} ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

if [ ${#xyz_pattern} -eq 0 ]
then
    # find the xyz roughly continuous
    myutils info_from_opt $pep # it looks for <pep>-opt.log and <pep>-stretched00.pdb
    xyz_pattern='-forces'
fi

cd $ds_dir
ds_path=$(pwd)
mapfile -t peptides < <(find . -maxdepth 1 -mindepth 1 -type d | sort)

for pep in ${peptides[@]}
do
  cd $pep
  if [ -d $subdir ]
  then
    cd $subdir
  fi
  # creates dofs. the outputs are *"xyz_pattern"*-dof.dat
  myutils extr_dofs -f "${pep#*/}"

  # send all selected files called <xyz_pattern>* to subset directory
  myutils reduce_structs "."
  cd subset
  myutils rearange_files
  cd $ds_path
done
