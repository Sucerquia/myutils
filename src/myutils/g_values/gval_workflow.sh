#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 16
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


#source $mine/sw/orca/setup_orca.sh

print_help() {
echo "
This tool runs all the necessary steps to get the g-values. It works with some
flags:
  
  -b  Use this flag to AVOID computation of frequencies.
  -c  <charge=0> charge of the system.
  -d  <directory='.'> directory containing at least model.xyz
  -e  Use this flag to AVOID computation of gvalues. Flags -f and -g are
      unuseful when this flag is activated.
  -f  <hyperfine=''> indexes of atoms to be included in hyperfine corrections.
      If the argument of this flag starts with g, the guess util
      'myutils HFC_relevantA' is used to find the H, O, and N atoms around the
      indexes atoms with the indexes given in the arguments. If you use the
      guess util, be sure to add the depth of the neighborhood with the letter
      d. E.g. 'g4,5d3' finds the relevant atoms in the 3-depth-neigborhood of
      the atoms 4 and 5.
  -g  <reference_mol=''> guess the location of the radical assuming an
      abstraction process. Give the path of the file of the molecule before
      the abstraction.
  -m  <multiplicity=2> multiplicity of the system.
  -o  Use this flag to AVOID the optimization step. opt.xyz must exist, then.
  -p  <processors=16> number of processors used in the orca calculations.

  -h   prints this message.
"
exit 0
}


charge=0
directory='.'
mult=2
optimization='true'
epr='true'
processors=16
hyperfine=''
bdes='true'

while getopts 'bc:d:ef:g:m:op:vh' flag;
do
  case "${flag}" in
    b) bdes='false' ;;
    c) charge=${OPTARG} ;;
    d) directory=${OPTARG} ;;
    e) epr='false' ;;
    f) hyperfine=${OPTARG} ;;
    g) guess_rad_loc=${OPTARG} ;;
    m) mult=${OPTARG} ;;
    o) optimization='false' ;;
    p) processors=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" Gvals $verbose

cd $directory

if $optimization
then
  verbose Optimization
  cat << EOF > opt.inp
! B3LYP EPR-II OPT
%pal nprocs $processors end
*XYZFile $charge $mult model.xyz
EOF
  $orca opt.inp  > opt.out
else
  [ -f opt.xyz ] || $orca opt.inp  > opt.out
fi

if $epr
then
  rad_loc=""
  if [[ "$guess_rad_loc" != "" ]]
  then
    location=$(myutils rad_loc $guess_rad_loc opt.xyz)
    rad_loc="$location,"
  fi

  if [[ ${hyperfine: 0: 1} == 'g' ]]
  then
    nog=${hyperfine: 1}           # remove g
    radicals="$rad_loc${nog%d*}"  # list of radicals
    depth=${nog#*d}               # depth
    tmp_var=$(myutils iHFC_fromxyz ${prior_name}_opt.xyz "[$radicals]" \
      "$depth")
    mapfile -t hyperfine < <(echo $tmp_var |  grep -oP '\[\K[^\]]+')
  fi

  verbose g-values
  cat <<EOF > epr_info.inp
! B3LYP EPR-II AUTOAUX
%pal nprocs $processors end
*XYZFile $charge $mult opt.xyz
%EPRNMR
        GTENSOR   TRUE
        ORI       GIAO
END
EOF
  # add necleous for hyperfine corrections
  if [[ $hyperfine != '' ]]
  then
    sed -i "/GTENSOR   TRUE/a\ \ \ \ \ \ \  NUCLEI\ \ \ \ = $hyperfine {SHIFT, AISO, ADIP, AORB}" epr_info.inp
  fi
  $orca epr_info.inp  > epr_info.out
fi

if $bdes
then
  verbose BDES
  myutils orca_epr_inp '' '' $processors $mult $charge 'opt.xyz' "$hyperfine" > epr_info.inp
  cat << EOF > freq.inp
! M062X def2-TZVP OPT FREQ
%pal nprocs $processors end
*XYZFile $charge $mult opt.xyz
EOF
  $orca freq.inp  > freq.out
fi
finish
