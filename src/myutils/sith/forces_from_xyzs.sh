#!/bin/bash

for pep in $(cat $1)
do
  echo ${pep%*/}
  cd $pep
  index1=$( grep ACE "${pep%*/}-stretched00.pdb" | grep CH3 | awk '{print $2}' )
  index2=$( grep NME "${pep%*/}-stretched00.pdb" | grep CH3 | awk '{print $2}' )

  cd forces
  for file in *forces*.xyz
  do
    tail -n +3 $file > tmp.xyz
    newzmat -ixyz -ozmat -rebuildzmat tmp ${file%.xyz}
    sed -i "s/-- No Title Specified --/Computation of forces/g" ${file%.xyz}.com
    sed -i "s/\# HF\/6-31G\* Test/%chk=${file%.xyz}\n%NProcShared=8\n#P bmk\/6-31+g opt(modredun,calcfc) force/g" ${file%*.xyz}.com
    echo -e "$index1 $index2 F" >>  ${file%.xyz}.com
  done
  rm tmp.xyz
  cd ../..
done
