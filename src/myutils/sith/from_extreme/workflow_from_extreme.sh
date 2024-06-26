#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


# ----- definition of functions starts ----------------------------------------
print_help() {
echo "
This tool creates the files to do the sith analysis by optimizing a molecule
that was just about to get a first rupture, then takes the intermedia steps and
find the internal forces. Consider the next options:

  -c  run in cascade. (modules are loaded)
  -p  <peptide>. directory or xyzfile of last conf.Chains of aminoacids to
      be evaluated. For example, \"./AAA/\" would optimize the last
      stretched a trialanine peptide (where last means after organizing
      alphabetically).
  -l  <number of amino acids in the peptide> It will be assumed that the
      xyz file starts with the letter code of the amino acids.

  -h  prints this message.
"
exit 0
}


resubmit () {
  sleep 23h 58m ; \
  if [[ "$(whoami)" == "hits_"* ]]
  then
    single_part="--partition=single"
  else
    single_part=""
  fi
  sbatch $single_part -J $SLURM_JOB_NAME "$( myutils workflow_from_extreme -path)" -p "$1" -c -r -s "$2"; \
  echo "new JOB submitted"
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
cascade='false'
n_processors=8
ref=''
restart=''
size=30
lenght=''

while getopts 'cl:p:rs:h' flag;
do
  case "${flag}" in
    c) cascade='true' ;;
    l) lenght=${OPTARG} ;;
    p) ref=${OPTARG} ;;
    r) restart='-r' ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" WF_FROM_EXTREME

if $cascade
then
    load_modules # ADD the parameters to resubmit
fi

# check that the variables are well set
if [ "${#lenght}" -eq 0 ]
then
  fail "You have to specify the lenght of the peptide usign the flag -l. For
    more information, use \"myutils workflow_from_extreme -h\""
fi

if [ "${#ref}" -eq 0 ]
then 
  fail "This code needs one reference. Please, define it using the flag -p.
    For more info, use \"myutils workflow_from_extreme -h\""
fi

# In case pep is a directory, it searches the last xyz in this dir.
if [ -d $ref ]
then
    cd $ref
    ref=${ref##*/}
    mapfile -t previous < <( find . -maxdepth 1 -type f -name "*$ref*.xyz" \
                                    -not -name "*bck*" | sort )
    ref=${previous[-1]}
fi

if [ ! -f $ref ]
then
    fail "$ref does not exist"
fi
# ---- BODY -------------------------------------------------------------------
xyz=${ref##*/}
name=${xyz:0:$lenght}

verbose "The first g09 process is an optimization starting from $ref"

# creates gaussian input that optimizes the structure
myutils change_distance "$xyz" "$name-optext" frozen_dofs.dat 0 0 \
  "scale_distance" || fail "Preparating the input of gaussian"
sed -i "1a %NProcShared=$n_processors" "$name-optext.com"
sed -i "3a opt(modredun,calcfc)" "$name-optext.com"

# run gaussian
verbose "Running optmization of stretching ${nameiplusone}"
g09 "$name-optext.com" "$name-optext.log" || \
  { if [ "$(grep -c "Atoms too close." \
         "$name-optext.com")" \
         -eq 1 ]; then fail "Atoms too close for ${nameiplusone}" ; \
    fi ; }

# check convergence from output
output=$(grep -i optimized "$name-optext.log" | \
            grep -c -i Non )

[ "$output" -ne 0 ] && fail "optimization didn't converged"

# TODO: add the smooter and the optimizator for the computation of forces.
finish "$name finish"
