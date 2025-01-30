#!/bin/bash

#SBATCH -N 1                   # number of nodes
#SBATCH -n 16
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --output=%x-%j.o
#SBATCH --error=%x-%j.e
#SBATCH --exclusive


source /hits/basement/mbm/sucerquia/sw/orca/setup_orca.sh
orca="/hits/basement/mbm/sucerquia/sw/orca/orca_5_0_4_linux_x86-64_openmpi411/orca"

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
basis='both'
charge=0

while getopts 'c:b:d:m:oh' flag;
do
    case "${flag}" in
      c) charge=${OPTARG} ;;
      b) basis=${OPTARG} ;;
      d) directory=${OPTARG} ;;
      m) mult=${OPTARG} ;;
      o) optimization='false' ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

reference=$(myutils gval_workflow -path)
reference=${reference%/basic_scripts*}

cd $directory
cp $reference/basic_scripts/*.inp .

for inp in *.inp
do
  sed -i "s/XYZFile 0 2/XYZFile $charge $mult/g" $inp
done

if [ "$basis" == "both" ]
then
  basis=( "EPRII" "cc-pVDZ" )
else
  basis=( $basis )
fi


for base in ${basis[@]}
do
  if $optimization
  then
    $orca opt_$base.inp  > opt_$base.out
  else
    [ -f opt_$base.xyz ] || $orca opt_$base.inp > opt_$base.out
  fi

  # Next line changed to delete base calculation and leave only base_i. Before: methods=( "${base}_i" "$base" )
  methods=( "${base}_i" )

  for met in ${methods[@]}
  do
      sbatch -J ${SLURM_JOB_NAME}_${met} $reference/basic_scripts/submit.sh "$met".inp "$met".out
  done
done
