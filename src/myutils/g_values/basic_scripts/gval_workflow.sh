#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 16
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


source $mine/sw/orca/setup_orca.sh

print_help() {
echo "
This tool runs all the necessary steps to get the g-values. It works with some flags:

  -c    <charge=0> charge of the system.
  -d    <directory='./'> directory containing at least model.xyz
  -m    <multiplicity=2> multiplicity of the system.
  -o    Use this flag to avoid the optimization step. opt.xyz must exist, then.

  -h   prints this message.
"
exit 0
}

directory='.'
mult=2
optimization='true'
epr='true'
charge=0
processors=16
hyperfine=''
bdes='true'

while getopts 'c:d:ef:m:op:vh' flag;
do
  case "${flag}" in
    b) bdes='false' ;;
    c) charge=${OPTARG} ;;
    d) directory=${OPTARG} ;;
    e) epr='false' ;;
    f) hyperfine=${OPTARG} ;;
    m) mult=${OPTARG} ;;
    o) optimization='false' ;;
    p) processors=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

source "$(myutils basics -path)" Gvals $verbose

reference=$(myutils gval_workflow -path)
reference=${reference%/basic_scripts*}

cd $directory

if $optimization
then
  verbose Optimization
  cat << EOF > opt_EPRII.inp
! B3LYP EPR-II OPT
%pal nprocs $processors end
*XYZFile $charge $mult model.xyz
EOF
  $orca opt_EPRII.inp  > opt_EPRII.out
else
  [ -f opt_EPRII.xyz ] || $orca opt_EPRII.inp  > opt_EPRII.out
fi

if [[ ${hyperfine: 0: 1} == 'g' ]]
then
  nog=${hyperfine: 1}
  radicals="${nog%d*}"
  depth=${nog#*d}
  tmp_var=$(myutils iHFC_fromxyz opt_EPRII.xyz "$radicals" "$depth")
  hyperfine="$(echo $tmp_var |\
             sed -E "s/\[//g ; s/\]//g; s/^ *//g ; s/ *$//g ; s/ +/,/g")"
fi

if $epr
then
  verbose g-values
  cat <<EOF > EPRII_i.inp
! B3LYP EPR-II AUTOAUX
%pal nprocs $processors end
*XYZFile $charge $mult opt_EPRII.xyz
%EPRNMR
        GTENSOR   TRUE
        ORI       GIAO
END
EOF
  # add necleous for hyperfine corrections
  if [[ $hyperfine != '' ]]
  then
    sed -i "/GTENSOR   TRUE/a        NUCLEI    = $hyperfine {SHIFT, AISO, ADIP, AORB}" EPRII_i.inp
  fi
  $orca EPRII_i.inp  > EPRII_i.out
fi

if $bdes
then
  verbose BDES
  myutils orca_epr_inp '' '' $processors $mult $charge 'opt_EPRII.xyz' "$hyperfine" > EPRII_i.inp
  cat << EOF > freq.inp
! M062X def2-TZVP OPT FREQ
%pal nprocs $processors end
*XYZFile $charge $mult opt_EPRII.xyz
EOF
  $orca freq.inp  > freq.out
fi
finish
