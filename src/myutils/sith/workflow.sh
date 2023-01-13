#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH --job-name="workflow"       #job name
#SBATCH -t 24:00:00
#SBATCH --output=peptides_analysis-%j.o
#SBATCH --error=peptides_analysis-%j.e
#SBATCH --exclusive


# set up starts
print_help() {
echo "
This tool executes the sith analysis in a set of peptides defined as arguments.
You may use this script as TEMPLATE to submit a Job in cascade or to execute
locally. Consider the next options:
    
    -a   chains of aminoacids to ve evaluated. For example, \"AAA\" would 
         analyse a trialanine peptide.
    -c   run in cascade. (modules are loaded)
    -n   pepgen options.
    -r   restart. In this case, run from the directory of the precreated peptide.

    -h   prints this message.

NOTE 1: if you do not give any flag, the code will follow the whole workflow,
namely, until aplying the JEDI method with SITH.
"
exit 0
}

# Function that returns the error message and stops the run if something fails.
function fail {
    printf '%s\n' "$1" >&2 
    exit "${2-1}" 
}

# General variables
peptides=$(echo $@ | sed "s/-.*//")
cascade='false'
restart=''


while getopts 'a:cnrh' flag; 
do
    case "${flag}" in
      a) peptides=${OPTARG} ;;
      c) cascade='true' ;;
      n) pep_options=${OPTARG} ;;
      r) restart='-r' ;;

      h) print_help
    esac
done

# starting information
echo "
    ++++++++ WorkFlow_MSG: VERBOSE - execution information +++++++++++++++++++"
echo " * Date:"
date
echo "* Command:"
echo $0 $@

if $cascade
then
    echo "This JOB will be run in the Node:"
    echo $SLURM_JOB_NODELIST
    cd $SLURM_SUBMIT_DIR
    # check dependencies
    module purge
    source $HOME/.bashrc
    source /etc/profile.d/modules.sh
    source /hits/basement/mbm/sucerquia/exec/load_g09.sh
fi

echo "
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

if [ ${#peptides} -eq 0 ]
then 
    fail "
    ++++++++ WorkFlow_MSG: ERROR - This code needs one peptide. Please, define it
    using the flag -a. For more info, use \"../workflow -h\" ++++++++"
fi

# set up finished

for pep in $peptides
do
    # Creation of the peptide directory and moving inside.

    # ---- firstly, backup previous directories with the same name
    if [[ $restart != '-r' ]]
    then
        pepgen -h &> /dev/null || fail "
    ++++++++ SYS_PULL_MSG: ERROR - This code needs pepgen ++++++++"
        bck=$pep-bck_1
        if [ -d $pep ]
        then 
            bck_i=2
            while [ -d $bck ]
            do
                bck=$pep-bck_$bck_i
                bck_i=$(( $bck_i + 1 ))
            done
            echo "
    ++++++++ WorkFlow_MSG: WARNING - pep directory arealdy exist. This directory
    will be backed up in $bck. ++++++++"
            mv $pep $bck
        fi
    
        mkdir $pep
        cd $pep
        # Creation of peptide
        pepgen $pep tmp -s flat $pep_options || fail "
    ++++++++ SYS_PULL_MSG: ERROR - Creating peptide $pep ++++++++"
        mv tmp/pep.pdb ./$pep-stretched00.pdb
        rm -r tmp
    else
        echo "
    ++++++++ WorkFlow_MSG: WARNING - restarted. ++++++++"
    fi
    
    # Stretching
    myutils stretching -p $pep $restart || fail "
    ++++++++ WorkFlow_MSG: ERROR - Stretching of $pep failed. ++++++++"

    myutils sith_analysis ./stretching/$pep-stretched00.fchk \
           ./stretching/ || fail "
    ++++++++ WorkFlow_MSG: ERROR - Sith analysis of $pep failed. ++++++++"
done

exit 0
