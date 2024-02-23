#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH --job-name="workflow"       #job name
#SBATCH -t 24:00:00
#SBATCH --output=peptides_analysis-%j.o
#SBATCH --error=peptides_analysis-%j.e
#SBATCH --exclusive


# ----- definition of functions starts ----------------------------------------
source "$(myutils basics -path)" WORKFLOW

print_help() {
echo "
This tool executes the sith analysis in a set of peptides defined as arguments.
You can use this code to submit a Job in cascade or to execute it locally. 
Consider the next options:
    
    -b    <number of breakages> The simulation will run until get this number of
              ruptures.
    -c    run in cascade. (modules are loaded)
    -d    <reference document> file containing the existing peptides to avoid
              repetition.
    -e    <endo> or <exo> states for initial state of proline. Default <random>.
    -m    <method>. stretching method. To see the options, use
              'myutils change_distance -h'
    -n    <options>. pepgen options.
    -p    <peptide>. Chains of aminoacids to be evaluated. For example, \"AAA\"
              would analyse a trialanine peptide.
    -R   random pepeptide. Give the number of amino acids with this argument.
    -r   restart. In this case, run from the directory of the precreated
         peptide.
    -s   <size[A]> of the step that increases the distances. Default 0.2A

    -h   prints this message.
"
exit 0
}

resubmit () {
    sleep 23h 58m ; \
    sbatch "$( myutils workflow -path)" -p "$1" -c -r -s "$2" -b "$3" -s "$4" ; \
    echo "new JOB submitted"
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
breakages=1
cascade='false'
endoexo='random'
method=0
pep=''
restart=''
size=0.2
ref_doc='00-aminos.txt'

while getopts 'd:b:ce:m:n:p:rR:s:h' flag;
do
    case "${flag}" in
      b) breakages=${OPTARG} ;;
      c) cascade='true' ;;
      d) ref_doc=${OPTARG} ;;
      e) endoexo=${OPTARG} ;;
      m) method=${OPTARG} ;;
      n) pep_options=${OPTARG} ;;
      p) pep=${OPTARG} ;;
      r) restart='-r' ;;
      R) random=${OPTARG} ;;
      s) size=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# random peptide
if [ ! "${#random}" -eq 0 ]
then
    [ -f "$ref_doc" ] || fail "Non-recognized $ref_doc, check flag -d"
    pep=$( myutils gen_randpep "$random" ) || fail "Creating random peptide"
    while awk '!/^#/ {print $1}' "$ref_doc" | grep -q "$pep"
    do
        pep=$( myutils gen_randpep "$random" ) || fail "Creating
            random peptide"
    done
    echo "$pep" "   R" >> "$ref_doc"
    warning "The code created the random peptide $pep, the workflow will run
        with this peptide even if you also passed -p argument."
fi

if [ "${#pep}" -eq 0 ]
then 
    fail "This code needs one peptide. Please, define it using the flag -p or
        -R. For more info, use \"myutils workflow -h\""
fi

if $cascade
then
    resubmit "$pep" "$method" "$breakages" "$size" &
    echo " * This JOB will be run in the Node:"
    echo "$SLURM_JOB_NODELIST"
    cd "$SLURM_SUBMIT_DIR" || fail "moving to execution directory: $SLURM_SUBMIT_DIR"
    source "$HOME/.bashrc"
    # shellcheck disable=SC1091
    source /hits/basement/mbm/sucerquia/exec/load_g09.sh
    conda activate myutils
    module purge
    module use /hits/sw/its/doserbd/haswell/modules/all/GROMACS
    module load 2020.3-fosscuda-2019b
fi

ase -h &> /dev/null || fail "This code needs ASE"
command -V g09 &> /dev/null || fail "This code needs gaussian"
gmx -h &> /dev/null || fail "This code needs gmx"
myutils -h &> /dev/null || fail "This code needs myutils"
perl -E "say '+' x 80"

# ----- set up finishes -------------------------------------------------------

# ---- firstly, backup previous directories with the same name
if [[ "$restart" != "-r" ]]
then
    # check pepgen
    pepgen -h &> /dev/null || fail "This code needs pepgen"

    # create back up
    create_bck "$pep"

    # Creation of the peptide directory and moving inside.
    mkdir "$pep"
    cd "$pep" || fail "directory $pep does not exist"
    verbose "generating peptide"
    # Creation of peptide
    # shellcheck disable=SC2086
    pepgen "$pep" tmp -s flat $pep_options || fail "Creating peptide $pep"
    mv tmp/pep.pdb "./$pep-stretched00.pdb"
    verbose "protonize"
    myutils protonize "./$pep-stretched00.pdb" "./$pep-stretched00.pdb" | \
        fail "protonizing"
    verbose "define proline state"
    myutils proline_mod -f "$pep-stretched00.pdb" -s "$endoexo" || \
        fail "Proline estates configuration"
    mv "$pep-stretched00modpro.pdb" "$pep-stretched00.pdb" 
    rm -r tmp
else
    # moving to the peptide directory
    cd "$pep" || fail "directory $pep does not exist"
    warning "$pep restarted"
fi

myutils stretching -b "$breakages" -p "$pep" "$restart" -m "$method" \
                   -s "$size" || fail "Stretching of $pep failed"

# Compute classical energies
verbose "computing classical energies."

myutils classical_energies
# compute forces
verbose "submitting comptutation of forces.";

# command -V nohub || fail "nohub does not exist"

module load slurm/20.11.7-1.hits
sbatch "$( myutils find_forces -path )" "tmp.out" &&
echo "computation of forces submitted"

finish "$pep finished"
exit 0
