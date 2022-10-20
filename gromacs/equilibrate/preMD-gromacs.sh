# This code tuns the steps to generate a simulation with GROMACS

# ==== Arguments ====
# 1 gromacs executer
# 2 Intial configuration as a pdb file
# 3 forcefield (see the list bellow)

gmx=$1
initial_conf=$2
forcefield=$3
nprocessors=$4
name=${initial_conf%.*}
eq_utils="/hits/basement/mbm/sucerquia/utils/equilibrate"

echo "ooo 0 - Removes water molecules in the protein" && \

ls && \
grep -v HOH $initial_conf > ${name}_clean.pdb && \

echo "ooo 1 - Creates the gmx files: ...top, ...gro, ...itp" && \
ls && \
echo "$forcefield\n" | $gmx pdb2gmx -f ${name}_clean.pdb -o ${name}_processed.gro -water spce &&\

echo "ooo 2 - Defines a box" && \
ls && \
$gmx editconf -f ${name}_processed.gro -o ${name}_newbox.gro -c -d 2.0 -bt cubic && \

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
mpirun -np $nprocessors $gmx mdrun -v -deffnm em  && \

echo "ooo 6 - Equilibration nvt"  && \
ls && \
$gmx grompp -f $eq_utils/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr  && \
mpirun -np $nprocessors $gmx mdrun -v -deffnm nvt  && \

echo "ooo 7 - Equilibration npt"  && \
ls && \
$gmx grompp -f $eq_utils/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr  && \
mpirun -np $nprocessors $gmx mdrun -deffnm npt  && \

echo "ooo 8 - finishes correctly" && \
ls

if [ 10 -eq 2 ] 
then
# THIS IS JUST A COMMENT WITH THE DIFFERENT FORCEFIELDS ALREADY IMPLEMENTED IN GROMACS,
# THE THIRD ARGUMENT WITH THIS IS EXECUTED MUST CORRESPOND TO ONE OF THESE NUMBERS.
 "1: AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24, 1999-2012, 2003)
 2: AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)
 3: AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29, 461-469, 1996)
 4: AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21, 1049-1074, 2000)
 5: AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65, 712-725, 2006)
 6: AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al., Proteins 78, 1950-58, 2010)
 7: AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)
 8: CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)
 9: GROMOS96 43a1 force field
10: GROMOS96 43a2 force field (improved alkane dihedrals)
11: GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)
12: GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)
13: GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)
14: GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856, DOI: 10.1007/s00249-011-0700-9)
15: OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)"
fi