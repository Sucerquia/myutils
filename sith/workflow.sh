#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH -p cascade.p
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
    
    -a   chains of aminoacids to ve evaluated. For example, \"AAA GG\" would 
         analyse a trialanine and a glycineglycine peptide.
    -c   run in cascade. (modules are loaded)
    -f   forces to be used in the pulling step in [kJ mol^-1 nm^-2]. 
         Default 300. The argument \"10 300 200\" define several forces.
    -g   gromacs binary. For example gmx or gmx_mpi. Default gmx.
    -p   follow the workflow until pulling classicaly and generating analysis.
    -o   follow the workflow until optimization of the stretched configuration.
    -s   follow the workflow until stretching by constrains.


    -h   prints this message.

NOTE 1: Please, add the flags after the peptides.

NOTE 2: if you do not give any flag, the code will follow the whole workflow,
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
utils="/hits/basement/mbm/sucerquia/utils/"
pulling='false'
opt='false'
stretch='false'
peptides=$(echo $@ | sed "s/-.*//")
cascade='false'
gmx='gmx'
forces="300" # [kJ mol^-1 nm^-2].


while getopts 'a:cg:opsf:h' flag; 
do
    case "${flag}" in
      a) peptides=${OPTARG} ;;
      c) cascade='true' ;;
      g) gmx=${OPTARG} ;;
      o) opt='true' ;;
      p) pulling='true' ;;
      s) stretch='true' ;;
      f) forces=${OPTARG} ;;

      h) print_help
    esac
done

echo "peptides $peptides"
if $cascade
then
    echo "This JOB will be run in the Node:"
    echo $SLURM_JOB_NODELIST
    cd $SLURM_SUBMIT_DIR
fi

# check dependencies
if [ ${#peptides} -eq 0 ]
then 
    fail "
    +++++ WorkFlow_MSG: ERROR - this code needs one peptide. Please, define at
    least one before flags +++++"
fi
# set up finished

for pep in $peptides
do
    # Creation of the peptide directory and moving inside.

    # ---- firstly, backup previous directories with the same name
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
    +++++ WorkFlow_MSG: WARNING - pep directory arealdy exist. This directory
    will be backed up in $bck. +++++"
        mv $pep $bck
    fi
    mkdir $pep
    cd $pep

    # Creation and pulling
    if $cascade
    then
    module purge
    module use /hits/sw/mcm/modules
    module load mcm/molbio/plumed/2.7.2-cascade
    export GMX_NO_MPI_BIN="/hits/fast/mcm/app/plumed/plumed-2.7.2/gromacs-2020.5/bin" 
    export OMP_NUM_THREADS=2
    $GMXBIN/gmx_mpi -h &> /dev/null && echo "++++ TEST Actually, it is working" || fail "
    +++++ TEEST: ERROR - this code needs gromacs ($gmx failed) +++++"
    $utils/gromacs/peptide_pulling.sh -p $pep -g $GMXBIN/gmx_mpi || fail "
    +++++ WorkFlow_MSG: ERROR - pulling of $pep failed. +++++"
    else
    $utils/gromacs/peptide_pulling.sh -p $pep -g $gmx || fail "
    +++++ WorkFlow_MSG: ERROR - pulling of $pep failed. +++++"
    fi
    
    if $pulling
    then 
        cd .. 
        echo "
    +++++ WorkFlow_MSG: WARNING - workflow stopped after pulling +++++"
        continue
    fi

    # Optimization
    if $cascade
    then
    module purge
    source $HOME/.bashrc
    source /etc/profile.d/modules.sh
    source /hits/basement/mbm/sucerquia/exec/load_g09.sh
    fi
    $utils/sith/optimization.sh -p $pep || fail "
    +++++ WorkFlow_MSG: ERROR - optimization of $pep failed. +++++"
    
    if $opt
    then 
        cd ..
        echo "
    +++++ WorkFlow_MSG: WARNING - workflow stopped after optimization +++++"
        continue
    fi

    # Stretching
    $utils/sith/stretching.sh -p $pep || fail "
    +++++ WorkFlow_MSG: ERROR - optimization of $pep failed. +++++"
    
    if $stretch
    then
        cd ..
        echo "
    +++++ WorkFlow_MSG: WARNING - workflow stopped after stretching +++++"
        continue
    fi

    # Going back to initial directory
    cd ..
done
