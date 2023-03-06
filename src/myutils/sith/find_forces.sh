#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --job-name="workflow"       #job name
#SBATCH -t 24:00:00
#SBATCH --output=compute_forces-%j.o
#SBATCH --error=compute_forces-%j.e
#SBATCH --exclusive

print_help() {
echo "
This tool computes the forces in all chk files and store them in a directory
called forces.

    -c    run in cascade
    -d    degrees of freedom. Use compute all the dofs in the opt file. By
          default, this code uses the DOFs of newzmat

    -h    prints this message.
"
exit 0
}

# Function that returns the error message and stops the run if something fails.
fail () {
    printf '%s\n' "$1" >&2
    exit "${2-1}"
}

compute_forces () {
    echo "++++ FORCES: construct Z-matrix for $1 ++++"
    newzmat -ichk -ozmat -rebuildzmat -bmodel $1 forces.com || fail "
    ++++ FORCES: error creating the matrix ++++"
    sed -i "s/#P bmk\/6-31+g opt(modredun,calcfc)/%NProcShared=8\n#P bmk\/6-31+g force/g" forces.com
    if [ -f extra_DOFs.dat ]
    then
        line=$( grep -n ".000" forces.com | tail -n 1 | cut -d ":" -f 1 )
        { head -n $line forces.com ; \
          cat extra_DOFs.dat ; \
          tail -n +$line forces.com ; } > tmp.com
        mv tmp.com forces.com
    fi
    g09 forces.com || fail "
    ++++ FORCES: computing forces ++++"
}

while getopts 'a:ch' flag; do
    case "${flag}" in
      a) all_directories=${OPTARG} ;;
      c) cascade='true' ;;

      h) print_help
    esac
done

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

echo "++++ FORCES: create directory ++++"
mkdir forces

chks=( $( ls *.chk ) )

echo "++++ FORCES: create extradofs.dat ++++"
compute_forces ${chks[0]}
$( myutils extract_forces ) || fail "
    ++++ FORCES: error extracting forces ++++"
myutils save_extradofs *00.fchk *00.fchk || fail "
    ++++ FORCES: error creating extradofs ++++"

for chkfile in ${chks[@]}
do
    compute_forces $chkfile
    name=$( echo $chkfile | sed "s/stretched/force/g" )
    echo "moving result to forces/${name%.*}.log"
    mv forces.log forces/${name%.*}.log
done

rm extra_DOFs.dat
rm forces.dat
mv forces.com forces/input_template.com
