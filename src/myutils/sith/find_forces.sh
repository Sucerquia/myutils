#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH --job-name="forces"       #job name
#SBATCH -t 24:00:00
#SBATCH --output=compute_forces-%j.o
#SBATCH --error=compute_forces-%j.e
#SBATCH --exclusive


# ----- definition of functions starts ----------------------------------------
source "$(myutils basics -path)" FIND_FORCES

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

compute_forces () {
    verbose "construct Z-matrix for $1"
    newzmat -ichk -ozmat -rebuildzmat -bmodel "$1" forces.com || fail "
    Error creating the matrix"
    sed -i "s/#P bmk\/6-31+g opt(modredun,calcfc)/%NProcShared=8\n#P bmk\/6-31+g force/g" forces.com
    verbose "executes g09 computation of forces for $1"
    g09 forces.com || fail "computing forces"
}

# ----- definition of functions finishes --------------------------------------

cascade='false'
directory='./'
pattern=''
while getopts 'd:cp:h' flag; do
    case "${flag}" in
      c) cascade='true' ;;
      d) directory=${OPTARG} ;;
      p) pattern=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if $cascade
then
    load_modules
fi

cd "$directory" || fail "moving to $directory"
verbose "Finding forces in the directory $( pwd )" 
verbose "Create forces directory and extracting forces"
create_bck forces
mkdir forces
mkdir bck
mv ./*-bck*.* bck

mapfile -t chks < <(ls "$pattern"*.chk)

for chkfile in "${chks[@]}"
do
    echo "$chkfile"
    compute_forces "$chkfile"
    name=${chkfile//stretched/force}
    verbose "Moving result to forces/${name%.*}.log"
    mv forces.log "forces/${name%.*}.log" || fail "moving results to forces
        directory."
done

mv forces.com forces/input_template.com || fail "moving template to forces
    directory"

myutils extract_forces

finish "finished"
exit 0
