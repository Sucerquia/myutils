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
print_help() {
echo "
This tool executes the sith analysis in a set of peptides defined as arguments.
You can use this code to submit a Job in cascade or to execute it locally. 
Consider the next options:
    
    -a   <peptide>. Chains of aminoacids to be evaluated. For example, \"AAA\" would 
         analyse a trialanine peptide.
    -c   run in cascade. (modules are loaded)
    -n   <options>. pepgen options.
    -R   random pepeptide. Give the number of amino acids with this argument.
    -r   restart. In this case, run from the directory of the precreated
         peptide.

    -h   prints this message.
"
exit 0
}

# Function that returns the error message and stops the run if something fails.
# Function that adjustes the text to 80 characters
adjust () {
    text=$( echo "++++++++ WORKFLOW_MSG: $@ " )
    addchar=$( expr 80 - ${#text} % 80 )
    text=$( echo $text $( perl -E "say '+' x $addchar" ))
    nlines=$( expr ${#text} / 80 )
    for (( w=0; w<=$nlines-1; w++ ))
    do
        echo ${text:$(( w * 79 )):79}
    done
    echo
}
# Function that returns the error message and stops the run if something fails.
fail () {
    adjust "ERROR" $1
    exit "${2-1}"
}

# prints some text adjusted to 80 characters per line, filling empty spaces
# with +
verbose () {
    adjust "VERBOSE" $1
}
warning () {
    adjust "WARNING" $1
}
resubmit () {
    sleep 23h 59m ; \
    sbatch $( myutils workflow ) -a $1 -c -r | echo ; \
    echo "new JOB submitted"
}
# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
pep=''
cascade='false'
restart=''
endoexo='random'

while getopts 'a:e:cn:rR:h' flag;
do
    case "${flag}" in
      a) pep=${OPTARG} ;;
      c) cascade='true' ;;
      e) endoexo=${OPTARG} ;;
      n) pep_options=${OPTARG} ;;
      r) restart='-r' ;;
      R) random=${OPTARG} ;;

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
    warning "The code will create a random peptide, even if you also passed -a
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
    fail "This code needs one peptide. Please, define it using the flag -a or
        -R. For more info, use \"myutils workflow -h\""
fi


if $cascade
then
    resubmit $pep &
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

# Stretching
$( myutils stretching ) -p $pep $restart || fail "Stretching of $pep failed"

# compute forces
$( myutils find_forces )

verbose "Workflow finished"
exit 0
