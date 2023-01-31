fail () {
    printf '%s\n' "$1" >&2
    exit "${2-1}"
}

counter=0

myutils all_xyz2pdb *-stretched00.pdb || fail "error creating pdbs"
echo "# counter bound angles dihedeal potential" > classical_energy.dat

for pdbfile in *.pdb
do
    echo $pdbfile
    for amino in PRO GLY
    do
        for leter in A B G D
        do
            for i in 2 3
            do
                sed -i "s/${i}H$leter  $amino/$(( $i -1 ))H$leter  $amino/g" $pdbfile
            done
        done
    done
    sed -i "s/ HA  GLY/1HA  GLY/g" $pdbfile
    echo -e "4\n 7\n" | gmx pdb2gmx -f $pdbfile -o minim.gro  || fail "error creating gro file"
    gmx editconf -f minim.gro -o minim_box.gro -c -d 5.0 -bt cubic &> /dev/null || fail "error creating box"
    mv minim_box.gro minim.gro || fail "error changing name"

    gmx grompp -f $( myutils minim ) -c minim.gro -p topol.top -o em.tpr &> /dev/null || fail "error creating em"

    gmx mdrun -v -deffnm em &> /dev/null || fail "error running em"

    echo -e "10 1 2 3 0\n" | gmx energy -f em.edr -o mini.xvg || fail "error computing energy"

    grep -v "@" mini.xvg | grep -v "#"| head -n 1 | \
    awk -v counter=$counter 'BEGIN{OFS="\t"}{print counter OFS $2 OFS $3 OFS $4 OFS $5}' \
    >> classical_energy.dat || fail "error saving energy"
    counter=$(( $counter + 1 ))
done

rm \#*
rm em.*
rm mdout.mdp
rm mini*
rm posre.itp
rm topol.top
