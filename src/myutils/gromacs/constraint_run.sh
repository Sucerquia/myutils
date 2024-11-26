#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Use this template to create your scripts with a standard structure

  -d  <variable> add the description of the variable.
  -c  usually used when submitted in a cluster to import modules.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
output='/dev/null'
gmx="gmx"
steps="100000"
verbose=''
cascade='false'
while getopts 'cd:g:l:s:vh' flag;
do
  case "${flag}" in
    d) distance=${OPTARG} ;;
    g) gmx=${OPTARG} ;;
    l) output=${OPTARG} ;;
    s) steps=${OPTARG} ;;
    c) cascade='true' ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" ConsDist $verbose

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# load modules
if $cascade
then
  echo "entering here"
  # Add flags to restart if necessary (and if you added a restart function before)
  load_modules 
fi

$gmx --version || fail "loading gmx"

# ---- BODY -------------------------------------------------------------------

distname=${distance//./}
mkdir E2Edist$distname
cd E2Edist$distname

cp "$( myutils constraint )" ./constraint.mdp && \
  sed -i "s/<ConsDist>/$distance/g" ./constraint.mdp && \
  sed -i "s/= 10000/= $steps/g" ./constraint.mdp || \
  fail "setting the file ./constraint.mdp"

$gmx grompp -f constraint.mdp \
            -c ../equilibrate/npt.gro \
            -t ../equilibrate/npt.cpt \
            -p ../equilibrate/pep_out.top \
            -n ../index.ndx \
            -maxwarn 5 \
            -o "md_0_$distname.tpr" > "$output" 2>&1 || \
  fail "grompp step of the constrain $distance"

verbose "MD run for distance $distance"

$gmx mdrun -deffnm "md_0_$distname" >> "$output" 2>&1 || \
  fail "Execution step of the constrain $distance"

rm -f \#*
myutils extract_distance -r 1,5 -a CH3,CH3 -t md_*.trr -o distance -g md_*.gro

finish "D=$distance finished"
