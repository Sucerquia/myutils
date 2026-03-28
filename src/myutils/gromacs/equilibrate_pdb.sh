#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --gres=gpu:1
#SBATCH --exclusive

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Use this template to create your scripts with a standard structure

  -f  <file> pdb file of the molecule to be pulled.
  -G  <name of the environment>Use it to use grappa in your simulations.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# source /etc/profile.d/z00_lmod.sh
# module load gmvolf/1.7.20
# ----- set up starts ---------------------------------------------------------
# General variables
cascade='false'
grappa='false'
verbose=''
gmx='gmx'
read_forces=''
steps=100000
while getopts 'f:g:G:cvh' flag;
do
  case "${flag}" in
    f) file=${OPTARG} ;;
    g) gmx=${OPTARG} ;;
    G) environment=${OPTARG} ;;
    c) cascade='true' ;;
    s) steps=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done


source "$(myutils basics -path)" equilibrate $verbose &&

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# load modules
if $cascade
then
  # Add flags to restart if necessary (and if you added a restart function before)
  load_modules 
fi

# ---- BODY -------------------------------------------------------------------
# ==== equilibrate
name=${file%.pdb}

mkdir equilibrate
cp $file equilibrate/pep.pdb
cd equilibrate

printf "4\n1\n " | $gmx pdb2gmx -f pep.pdb \
                                -o pep.gro \
                                -p pep_out.top -ignh \
  || fail "pdb2gmx failed"

if [[ ${#environment} -ne 0 ]]
then
  verbose "Activate grappa and construct topology"
  #CONDA_PREFIX="/hits/basement/mbm/sucerquia/conda"
  #eval "$($CONDA_PREFIX/bin/conda shell.bash hook)"
  #conda init
  #conda activate $environment 
  grappa_gmx -h > /dev/null || fail "grappa is not installed. Be sure it is installed in
    the selected environment."
  grappa_gmx -f pep_out.top -o topology_grappa.top -t grappa-1.3.0 \
    || fail "Error of grappa creating the topology"
  mv topology_grappa.top pep_out.top
fi

verbose "create simulation cell"
$gmx editconf -f pep.gro \
              -o pep_box.gro \
              -c -d 1.0 -bt dodecahedron \
  || fail "editconfig failed"

# Solvate
verbose solvate
$gmx solvate -cp pep_box.gro \
             -p pep_out.top \
             -o pep_solv.gro || fail "solvate failed"

# add ions
verbose "add ions"
cp "$( myutils ions )" ./ions.mdp
$gmx grompp -f ions.mdp \
            -c pep_solv.gro \
            -p pep_out.top \
            -o pep_out_genion.tpr \
  || fail "grompp ions failed"
echo "SOL" | $gmx genion -s pep_out_genion.tpr \
                         -p pep_out.top \
                         -o pep_out_ion.gro \
                         -conc 0.15 -neutral \
  || fail "genion failed"

# minimize
verbose "minimize"
cp "$( myutils minim )" ./minim.mdp
$gmx grompp -f minim.mdp \
            -c pep_out_ion.gro \
            -p pep_out.top \
            -o pep_out_min.tpr \
  || fail "grompp minimization failed"

$gmx mdrun -deffnm pep_out_min -v \
  || fail "minimization failed"

# nvt
verbose "nvt"
cp "$( myutils nvt )" ./nvt.mdp
$gmx grompp -f nvt.mdp \
            -c pep_out_min.gro \
            -p pep_out.top \
            -o nvt.tpr \
  || fail "grompp nvt failed"
$gmx mdrun -v -deffnm nvt \
  || fail "nvt failed"

# npt
verbose "npt"
cp "$( myutils npt )" ./npt.mdp
$gmx grompp -f npt.mdp \
            -c nvt.gro \
            -p pep_out.top \
            -o npt.tpr \
  || fail "grompp npt failed"
$gmx mdrun -v -deffnm npt \
  || fail "npt failed"

# going back to the original directory
cd ..

finish
