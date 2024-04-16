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
while getopts 'd:s:h' flag;
do
  case "${flag}" in
    d) ds_dir=${OPTARG} ;;
    s) subdir=${OPTARG} ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

original_path=$(pwd)
cd $ds_dir
ds_path=$(pwd)
mapfile -t peptides < <(find . -maxdepth 1 -mindepth 1 -type d | sort)

for pep in ${peptides[@]}
do
  echo $pep
  cd $pep
  if [ -d $subdir ]
  then
    cd $subdir
  fi
  # creates dofs. the outputs are ${<*pep*.xyz>%.xyz}-dof.dat
  myutils extr_dofs -f "${pep#*/}" > /dev/null

  # send all selected files called ${<*pep*.xyz>%.xyz}* to subset directory
  myutils reduce_structs "."
  echo -n "rename: "
  cd subset
  myutils rearange_files
  cd $ds_path
done

cd $original_path
