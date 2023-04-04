#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH --job-name="workflow"       #job name
#SBATCH -t 24:00:00
#SBATCH --output=compute_forces-%j.o
#SBATCH --error=compute_forces-%j.e
#SBATCH --exclusive

print_help() {
echo "
This tool computes the forces in all chk files and store them in a directory
called forces.

    -c    run in cascade.
    -d    directory containging the chk files of the stretching-optimization
          process. Default ./

    -h    prints this message.
"
exit 0
}

# Function that adjustes the text to 80 characters
adjust () {
    text=$( echo "++++++++ FORCES: $@ " )
    addchar=$( expr 80 - ${#text} % 80 )
    text=$( echo $text $( perl -E "say '+' x $addchar" ))
    nlines=$( expr ${#text} / 80 )
    w=0
    while [ $w -le $(( nlines - 1 )) ]
    do
        echo ${text:$(( w * 79 )):79}
        w=$(( w + 1 ))
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

compute_forces () {
    verbose "construct Z-matrix for $1"
    newzmat -ichk -ozmat -rebuildzmat -bmodel $1 forces.com || fail "
    Error creating the matrix"
    sed -i "s/#P bmk\/6-31+g opt(modredun,calcfc)/%NProcShared=8\n#P bmk\/6-31+g force/g" forces.com
    g09 forces.com || fail "computing forces"
}

# ----- definition of functions finishes --------------------------------------

directory='./'
while getopts 'd:ch' flag; do
    case "${flag}" in
      d) directory=${OPTARG} ;;
      c) cascade='true' ;;

      h) print_help
    esac
done
cd $directory

if $cascade
then
    echo " * This JOB will be run in the Node:"
    echo $SLURM_JOB_NODELIST
    cd $SLURM_SUBMIT_DIR
    # check dependencies
    module purge
    source $HOME/.bashrc
    source /etc/profile.d/modules.sh
    source /hits/basement/mbm/sucerquia/exec/load_g09.sh
    conda activate myutils
fi

verbose "Create forces directory"
mkdir forces
mkdir bck
mv *-bck*.* bck
mv *-a.* bck

chks=( $( ls *.chk ) )

for chkfile in ${chks[@]}
do
    compute_forces $chkfile
    name=$( echo $chkfile | sed "s/stretched/force/g" )
    verbose "Moving result to forces/${name%.*}.log"
    mv forces.log forces/${name%.*}.log
done

rm forces.dat
mv forces.com forces/input_template.com
