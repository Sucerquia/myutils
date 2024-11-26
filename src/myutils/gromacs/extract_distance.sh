#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 9
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Extract distance between two atoms from gromacs trajectory.

  -r  <res1,res2>, residues indexes of the atoms to compute the distance.
  -a  <a1, a2>, names of the atoms to compute the distance.
  -t  <traj file> path to the trajectory file to extract the distance.
  -o  <output> name of the output without extension (.dat)
  -g  <gro_file> gro file used to create the trajectory.

  -v  verbose.
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
r1=''
r2=''
a1=''
cascade='false'
verbose=''
stream='/dev/null'
while getopts 'r:a:t:o:g:cs:vh' flag;
do
  case "${flag}" in
    r) res=${OPTARG} ;;
    a) atoms=${OPTARG} ;;
    t) traj=${OPTARG} ;;
    o) out=${OPTARG} ;;
    g) gro=${OPTARG} ;;
    c) cascade='true' ;;
    s) stream=${OPTARG} ;;

    v) verbose='-v' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" ExtDistance $verbose

if [ "$verbose" = '-v' ]
then
  stream='/dev/tty'
fi

# starting information
verbose " Command: " "$0" "$@"

# load modules
if $cascade
then
  # Add flags to restart if necessary (and if you added a restart function before)
  load_modules 
fi

r1=$( echo $res | cut -d ',' -f 1 )
r2=$( echo $res | cut -d ',' -f 2 )
if [ "$r1" == '' ]
then
  fail "define residues using the flag -r res1_index,res2_index. E.g. -r 1,2"
fi

a1=$( echo $atoms | cut -d ',' -f 1 )
a2=$( echo $atoms | cut -d ',' -f 2 )
if [ "$a1" == "" ]
then
  fail "define atoms using the flag -a atom1_name,atom2_name. E.g. -a Ca,C"
fi

# ---- BODY -------------------------------------------------------------------
verbose "Create index file (${a1}${r1}_${a2}${r2}.ndx)"

echo -e "ri $r1 & a $a1 \n ri $r2 & a $a2 \n \"r_${r1}_&_${a1}\" | \"r_${r2}_&_${a2}\" \n q \n" \
  | gmx make_ndx -f $gro -o ${a1}${r1}_${a2}${r2}.ndx > $stream 2>&1 || fail "make index"

name="${r1}${a1}_${r2}${a2}"
sed -i "s/r_${r1}_&_${a1}_r_${r2}_&_${a2}/$name/g" ${a1}${r1}_${a2}${r2}.ndx || \
  fail "renaming group in index file"

verbose "compute distance"

echo -e "\"$name\"\n" | \
  gmx distance -f $traj -s $gro -n ${a1}${r1}_${a2}${r2}.ndx -oall "$out.xvg" > \
  $stream 2>&1 || fail "computation of distance"
grep -v "^#" "$out.xvg" > tmp.dat
grep -v "^@" tmp.dat > "$out.dat"
sed -i "1i # time[ns] distance[nm]" "$out.dat"

finish
