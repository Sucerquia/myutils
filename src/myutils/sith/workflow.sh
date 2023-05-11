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
source $(myutils basics) STRETCHING_MSG

print_help() {
echo "
This tool executes the sith analysis in a set of peptides defined as arguments.
You can use this code to submit a Job in cascade or to execute it locally. 
Consider the next options:
    
    -b    <number of breakages> The simulation will run until get this number of
              ruptures.
    -c    run in cascade. (modules are loaded)
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
    sleep 23h 59m ; \
    sbatch $( myutils workflow ) -a $1 -c -r -s $2 | echo ; \
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

while getopts 'b:ce:m:n:p:rR:s:h' flag;
do
    case "${flag}" in
      b) breakages=${OPTARG} ;;
      c) cascade='true' ;;
      e) endoexo=${OPTARG} ;;
      m) method=${OPTARG} ;;
      n) pep_options=${OPTARG} ;;
      p) pep=${OPTARG} ;;
      r) restart='-r' ;;
      R) random=${OPTARG} ;;
      s) size=${OPTARG} ;;

      h) print_help
    esac
done

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo $0 $@

# random peptide
if [ ! ${#random} -eq 0 ]
then
    warning "The code will create a random peptide, even if you also passed -p
        argument"
    pep=$( myutils gen_randpep $random ) || fail "Creating
        random peptide"
    while [ -d $pep ]
    do
        pep=$( myutils gen_randpep $random ) || fail "Creating
            random peptide"
    done
fi

if [ ${#pep} -eq 0 ]
then 
    fail "This code needs one peptide. Please, define it using the flag -p or
        -R. For more info, use \"myutils workflow -h\""
fi


if $cascade
then
    resubmit $pep $stretching_type &
    echo " * This JOB will be run in the Node:"
    echo $SLURM_JOB_NODELIST
    cd $SLURM_SUBMIT_DIR
    # check dependencies
    source $HOME/.bashrc
    source /hits/basement/mbm/sucerquia/exec/load_g09.sh
    conda activate myutils
    module purge
    module use /hits/sw/its/doserbd/haswell/modules/all/GROMACS
    module load 2020.3-fosscuda-2019b

    ase -h &> /dev/null || fail "This code needs ASE"
    command -V g09 &> /dev/null || fail "This code needs gaussian"
    gmx -h &> /dev/null || fail "This code needs gmx"
    myutils -h &> /dev/null || fail "This code needs myutils"
fi
perl -E "say '+' x 80"

# ----- set up finishes -------------------------------------------------------

# ---- firstly, backup previous directories with the same name
if [[ $restart != '-r' ]]
then
    # check pepgen
    pepgen -h &> /dev/null || fail "This code needs pepgen"

    # create back up
    bck=$pep-bck_1
    if [ -d $pep ]
    then 
        bck_i=2
        while [ -d $bck ]
        do
            bck=$pep-bck_$bck_i
            bck_i=$(( $bck_i + 1 ))
        done
        warning "$pep directory already exist. This directory will be
            backed up in $bck"
        mv $pep $bck
    fi
    # Creation of the peptide directory and moving inside.
    mkdir $pep
    cd $pep
    # Creation of peptide
    pepgen $pep tmp -s flat $pep_options || fail "Creating peptide $pep"
    mv tmp/pep.pdb ./$pep-stretched00.pdb
    $(myutils proline_mod) $pep-stretched00.pdb $endoexo

    rm -r tmp
else
    # moving to the peptide directory
    cd $pep
    warning "restarted"
fi

$( myutils stretching ) -b $breakages -p $pep $restart -m $method \ 
    -s $size || fail "Stretching of $pep failed"

# Compute classical energies
verbose "computing classical energies."

$( myutils classical_energies )
# compute forces
verbose "submitting comptutation of forces."
pwd
sbatch $( myutils find_forces ) -c

verbose "Workflow finished"

finish
exit 0
