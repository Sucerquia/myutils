#!/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics) peptide_pulling

print_help() {
echo "
This tool creates the trajectory of a given peptide pulled by an external force.
Consider the next options:
   
    -a    properties you want to analyse. For example \"-d -r\". Default \"-d -l\".
          For more information, check: myutils analysis -h
    -f    forces to stretch the peptide in [kJ mol^-1 nm^-1]. Default 200
    -g    gromacs binary. For example gmx or gmx_mpi. Default gmx.
    -o    pepgen flags
    -p    peptide.

    -h    prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
gmx='gmx'
forces=()
analysis="-d -L"
steps="10000"
pep_options=""

while getopts 'a:f:g:op:s:h' flag; do
    case "${flag}" in
      a) analysis=${OPTARG} ;;
      f) forces=${OPTARG} ;;
      g) gmx=${OPTARG} ;;
      o) pep_options=${OPTARG} ;;
      p) pep=${OPTARG} ;;
      s) steps=${OPTARG} ;;

      h) print_help
    esac
done

# check dependencies
pepgen -h &> /dev/null || fail "This code needs pepgen"
$gmx -h &> /dev/null || fail "This code needs gromacs ($gmx failed)"
# ----- set up finishes -------------------------------------------------------

# create peptide
verbose "Creation and equilibration of $pep starts"
pepgen $pep equilibrate -gmx $gmx $pep_options -e || fail "Creating peptide $pep"

# pulling
verbose "Pulling of $pep starts"
if [ ${#forces[@]} -eq 0 ]
then
    warning "You didn't specify which forces you want to use. Then the pulling
        will be done using the value by default: 200 [kJ mol^-1 nm^-1]"
    forces="200"
fi

for force in $forces
do
    forcename=$(printf "%04d" $force)
    verbose "Force $force acting $pep  starts"
    $( myutils pulling ) -g $gmx -f $force -s $steps || fail "Pulling $pep with
        $force failed"
    
    create_bck force$forcename
    mkdir force$forcename && \
        mv md_0_* force$forcename && \
        mv *.ndx force$forcename && \
        mv *.mdp force$forcename || fail "Moving gromacs files to the ditectory
           force$forcename"

    verbose "Analysis of $pep and $force starts"
    cd force$forcename
    file=$( ls md_*.gro)
    name=${file%%.*}
    $( myutils analysis ) -f $name -g $gmx $analysis || \
        fail "Equilibration Analysis."
    cd ..
done

verbose "$pep pulling finishes"
finish
exit 0
