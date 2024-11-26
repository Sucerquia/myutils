#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --gres=gpu:1
#SBATCH --exclusive

# This code submit an equilibration and run a pulling simulation using
# grappa and amber99.
# TODO: I have to add this part to the workflow such that I don't
# have to run it afterwards.

file=$1
cascade=$2

source "$(myutils basics -path)" PULLING -v

if [ ${#cascade} -eq 0 ]
then
  cascade=''
else
  source /etc/profile.d/z00_lmod.sh
  load_modules
  cascade='sbatch'
fi
  

stretching() {
  force_int=$1
  if [ -d "force${force_int%.*}" ] && [[ $(ls "force${force_int%.*}/" | wc -l) -ne 0 ]]
  then
    echo "force${force_int%.*} exists"
  else
    echo running "force${force_int%.*}"
    rm -rf "force${force_int%.*}"
    $(myutils pulling -path) -f $force_int -v -l "/dev/tty" \
      || fail "pulling from $(pwd)"
    mkdir force${force_int%.*}
    echo "mv md_0_${force_int%.*}* force${force_int%.*}"
    mv md_0_${force_int%.*}* force${force_int%.*}
  fi
}

mapfile -t aas < <( awk '!/^#/ {print $1}' $file )
mapfile -t ffs < <( awk '!/^#/ {print $2}' $file )

for (( i=0; i<${#ass[@]}; i++ ))
do
  amino=${aas[i]}
  force=${ffs[i]}
  cd $amino

  # Grappa
  echo -e "$i $amino GRAPPA"
  mkdir grappaforced
  cp $amino-stretched00.pdb grappaforced || fail "copying pdb"
  cd grappaforced
  if [ ! -f "equilibrate/npt.trr" ]
  then
    $(myutils equilibrate_pdb -path) -f $amino-stretched00.pdb -v -G grappa \
      || fail "equilibrating $amino with grappa"
  fi

  stretching "$force" || fail "pulling from $(pwd) $force"
  stretching 3000 || fail "pulling from $(pwd) 3000"
  cd ..

  # amber99
  echo -e "\n\n\n $i $amino amber99 \n\n\n"
  mkdir amber99forced
  cp $amino-stretched00.pdb amber99forced 
  cd amber99forced
  if [ ! -f "equilibrate/npt.trr" ]
  then
    $(myutils equilibrate_pdb -path) -f $amino-stretched00.pdb -v \
      || fail "equilibrating $amino with amber"
  fi
  stretching "$force" || fail "pulling from $(pwd) $force"
  stretching 3000 || fail "pulling from $(pwd) 3000"

  cd ../../
done
