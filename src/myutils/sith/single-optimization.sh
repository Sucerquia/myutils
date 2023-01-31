#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 8
#SBATCH --job-name="optimization"       #job name
#SBATCH -t 24:00:00
#SBATCH --output=optimization-%j.o
#SBATCH --error=optimization-%j.e
#SBATCH --exclusive

if [ $1 == '-h' ]
then
    echo "
    This code runs one optimization using gaussian. You have to create the
    input file and give it as first argument when run this code."

echo "This JOB will be run in the Node:"
echo $SLURM_JOB_NODELIST
cd $SLURM_SUBMIT_DIR
# check dependencies

module purge
source $HOME/.bashrc
source /etc/profile.d/modules.sh
source /hits/basement/mbm/sucerquia/exec/load_g09.sh

g09 $1.com $1.log

