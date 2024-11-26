#!/bin/bash

# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
Creates the com files from the xyz structures extracted from a g09 log file and
submit the corresponding jobs to compute the forces.

  -l  <log_file> optimization g09 logfile.
  -n  <name> standard name. Usually pep name.

  -h  prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
while getopts 'l:n:h' flag;
do
  case "${flag}" in
    l) logfile=${OPTARG} ;;
    n) name=${OPTARG} ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source $(myutils basics -path) "after_opt"
# ---- BODY -------------------------------------------------------------------

# ==== Reduce number of structures with reduced changes of DOFs
verbose "Create continuous structutes path."
echo "Removes high energy"
# The output are the xyz files without peak energies, output name-forces<n>.xyz
myutils info_from_opt $logfile ${name}-stretched00.pdb ${name}-forces > \
  /dev/null || fail "extracting xyz files from log file"
# Extract the dofs from the created xyzs. out; <name>-forces-dofs.dat
myutils extr_dofs -f ${name}-forces > /dev/null || \
  fail "extracting dofs from xyzs"
# reduce irrelevant changes, store the new subset in a dir called subset
myutils reduce_structs "." ${name}-forces > /dev/null || \
  fail "reducing structures"

# ==== Create com g09 files
# Create com file template
myutils forces_from_xyzs -d . -n ${name}-forces000 -p ${name}-stretched00.pdb \
  > /dev/null|| fail "creating com files using forces_from_xyz"
# clean files: only leaves the template
mv ${name}-forces000.com template.com
sed -i "s/opt(modredun,calcfc) //g" template.com
echo "" >> template.com
rm *forces*

# import xyz files of the subset
mv subset/* .
rm -r subset

# Create .com files
verbose "Create com files."
str_index=0
for file in ${name}-forces*.dat
do
  struct_name=${file%.dat}
  echo $struct_name
  myutils find_blocks -f template.com -e "Variables:" -o tmp > /dev/null
  mv tmp_001.out  $struct_name.com
  echo "     Variables:" >> $struct_name.com
  cat $file >> $struct_name.com
  echo "" >> $struct_name.com
  myutils find_blocks -s "\^\$" -e "\^\$" -f template.com -o tmp > /dev/null
  cat tmp_001.out >> $struct_name.com
  sbatch -J ${name}_f-$str_index \
         $(myutils single_g09 -path) -f $struct_name \
	                                   -c || fail "submitting forces Job"
  str_index=$(( 10#$str_index + 1 ))
done

rm tmp_001.out
rm template.com
rm *.dat
rm *.pdb

finish
