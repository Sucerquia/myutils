# This code tuns the steps to generate a simulation with GROMACS

# ==== Arguments ====
# 1 gromacs executer
# 2 Intial configuration as a pdb file
# 3 forcefield (see the list bellow)

gmx=$1
initial_conf=$2
forcefield=$3
name=${initial_conf%.*}
eq_utils="/hits/basement/mbm/sucerquia/utils/equilibrate"

echo "ooo 0 - Removes water molecules in the protein" && \

ls && \
grep -v HOH $initial_conf > ${name}_clean.pdb && \

echo "ooo 1 - Creates the $gmx files: ...top, ...gro, ...itp" && \
ls && \
echo "$forcefield\n" | $gmx pdb2gmx -f ${name}_clean.pdb -o ${name}_processed.gro -water spce &&\

echo "ooo 2 - Defines a box" && \
ls && \
$gmx editconf -f ${name}_processed.gro -o ${name}_newbox.gro -c -d 1.0 -bt cubic && \

echo  "ooo 3 - Fills the box with water molecules" && \
ls && \
$gmx solvate -cp ${name}_newbox.gro -cs spc216.gro -o ${name}_solv.gro -p topol.top  && \

echo "ooo 4 - Generates a file with ions"  && \
ls && \
$gmx grompp -f $eq_utils/ions.mdp -c ${name}_solv.gro -p topol.top -o ions.tpr  && \
echo "13\n" | $gmx genion -s ions.tpr -o ${name}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral  && \

echo "ooo 5 - Minimization to start with a stable configuration"  && \
ls && \
$gmx grompp -f $eq_utils/minim.mdp -c ${name}_solv_ions.gro -p topol.top -o em.tpr  && \
mpirun -np 16 $gmx mdrun -v -deffnm em  && \

echo "ooo 6 - Equilibration nvt"  && \
ls && \
$gmx grompp -f $eq_utils/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr  && \
mpirun -np 16 $gmx mdrun -v -deffnm nvt  && \

echo "ooo 7 - Equilibration npt"  && \
ls && \
$gmx grompp -f $eq_utils/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr  && \
mpirun -np 16 $gmx mdrun -deffnm npt  && \

echo "ooo 8 - finishes correctly"
ls && \
