#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 16
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


source /hits/basement/mbm/sucerquia/sw/orca/setup_orca.sh

/hits/basement/mbm/sucerquia/sw/orca/orca_5_0_4_linux_x86-64_openmpi411/orca $1 > $2
