#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH --job-name="optimization"       #job name
#SBATCH -t 24:00:00
#SBATCH --output=optimization-%j.o
#SBATCH --error=optimization-%j.e
#SBATCH --exclusive


print_help() {
echo "
This code runs one optimization using gaussian in one of the clusters. You
have to create the input file and give it (without .com extension) as first
argument when run this code.

    -c    run in cascade.

    -h    prints this message.
"
exit 0
}

cascade='false'
while getopts 'd:cp:h' flag; do
    case "${flag}" in
      c) cascade='true' ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if $cascade
then
    echo "This JOB will be run in the Node:"
    echo "$SLURM_JOB_NODELIST"
    cd "$SLURM_SUBMIT_DIR" || fail "moving to the directory from where the job was 
        submitted $SLURM_SUBMIT_DIR"

    module purge

    source "$HOME/.bashrc"
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
    # shellcheck disable=SC1091
    source /hits/basement/mbm/sucerquia/exec/load_g09.sh
fi

g09 "$1.com" "$1.log"
