#!/bin/bash

source "$(myutils basics -path)" PrepareAndSubmit

for pep in */
do
    verbose ${pep%*/}
    cd $pep
    index1=$( grep ACE "${pep%*/}-stretched00.pdb" | grep CH3 | awk '{print $2}' )
    index2=$( grep NME "${pep%*/}-stretched00.pdb" | grep CH3 | awk '{print $2}' )

    # check that the indexes were read properly:
    [[ "$index1" -eq 0 && "$index2" -eq 0 ]] && fail "Not recognized indexes"
    [[ "$index1" -eq 1 && "$index2" -eq 1 ]] && fail "Not recognized indexes"
    echo "$pep freezing $index1 and $index2"

    cd forces
    for file in *.xyz
    do
        myutils change_distance \
                $file ${file%.xyz}-opt \
                "nofile" 0 0 "scale_distance" \
                || fail "Preparating g09 input"
        sed -i '$d' ${file%.xyz}-opt.com
        echo "$index1 $index2 F" >> ${file%.xyz}-opt.com
        sed -i "1a %NProcShared=8" "${file%.xyz}-opt.com"
        sed -i "3a opt(modredun,calcfc)" "${file%.xyz}-opt.com"
        sbatch --job-name="${file:0:6}_opt" \
               --output="${file:0:6}_opt.o" \
               --error="${file:0:6}_opt.e" \
               $(myutils opt_and_forces -path) -f ${file%.xyz}-opt -c
    done
    cd ../..
done


