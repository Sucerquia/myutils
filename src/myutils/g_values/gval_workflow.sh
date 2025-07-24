#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 16
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e


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

bdes='true'
charge=0
directory='.'
mult=2
optimization='true'
epr='true'
processors=16
hyperfine=''
prior_name='model'

while getopts 'bc:d:ef:m:n:op:r:svh' flag;
do
  case "${flag}" in
    b) bdes='false' ;;
    c) charge=${OPTARG} ;;
    d) directory=${OPTARG} ;;
    e) epr='false' ;;
    f) hyperfine=${OPTARG} ;;
    m) mult=${OPTARG} ;;
    n) prior_name=${OPTARG} ;;
    o) optimization='false' ;;
    p) processors=${OPTARG} ;;
    r) reference_mol=${OPTARG} ;;
    s) sweep='true' ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" Gvals $verbose
load_modules

cd $directory

# ==== optimization ===========================================================
if $optimization
then
  verbose Optimization
  cat << EOF > ${prior_name}_opt.inp
! B3LYP EPR-II OPT
%pal nprocs $processors end
*XYZFile $charge $mult ${prior_name}.xyz
EOF
  $orca ${prior_name}_opt.inp  > ${prior_name}_opt.out
else
  [ -f ${prior_name}_opt.xyz ] || $orca ${prior_name}_opt.inp  > \
    ${prior_name}_opt.out
fi

if grep -q "ORCA TERMINATED NORMALLY" ${prior_name}_opt.out && $sweep
then
  rm ${xyz_file%.xyz}_opt.densities
  rm ${xyz_file%.xyz}_opt.engrad
  rm ${xyz_file%.xyz}_opt.gbw
  rm ${xyz_file%.xyz}_opt.inp
  rm ${xyz_file%.xyz}_opt.opt
  rm ${xyz_file%.xyz}_opt_property.txt
  rm ${xyz_file%.xyz}_opt_trj.xyz
  rm ${xyz_file%.xyz}_opt.gori.xyz
fi

# ==== epr ====================================================================
if $epr
then
  rad_loc=""
  if [[ "$reference_mol" != "" ]]
  then
    location=$(myutils rad_loc $reference_mol ${prior_name}_opt.xyz)
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
  cat <<EOF > ${prior_name}_epr.inp
! B3LYP EPR-II AUTOAUX
%pal nprocs $processors end
*XYZFile $charge $mult ${prior_name}_opt.xyz
%EPRNMR
        GTENSOR   TRUE
        ORI       GIAO
END
EOF
  # add necleous for hyperfine corrections
  if [[ $hyperfine != '' ]]
  then
    verbose hyperfine
    for sublist in "${hyperfine[@]}"
    do
      sed -i "/GTENSOR   TRUE/a\ \ \ \ \ \ \  NUCLEI\ \ \ \ = \
        $sublist {SHIFT, AISO, ADIP, AORB}" ${prior_name}_epr.inp
    done
  fi
  $orca ${prior_name}_epr.inp  > ${prior_name}_epr.out
fi

if grep -q "ORCA TERMINATED NORMALLY" ${prior_name}_epr.out && $sweep
then
  rm ${xyz_file%.xyz}_epr.densities
  rm ${xyz_file%.xyz}_epr.engrad
  rm ${xyz_file%.xyz}_epr.gbw
  rm ${xyz_file%.xyz}_epr.inp
  rm ${xyz_file%.xyz}_epr.opt
  rm ${xyz_file%.xyz}_epr_property.txt
  rm ${xyz_file%.xyz}_epr_trj.xyz
  rm ${xyz_file%.xyz}_epr.gori.xyz
fi

# ==== BDEs ===================================================================
if $bdes
then
  verbose BDES
  cat << EOF > ${prior_name}_freq.inp
! M062X def2-TZVP OPT FREQ
%pal nprocs $processors end
*XYZFile $charge $mult ${prior_name}_opt.xyz
EOF
  $orca ${prior_name}_freq.inp  > ${prior_name}_freq.out
fi

finish
