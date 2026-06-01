#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Create all the files after complete the computation of the gvalues with Orca.

  -c  <molecule candidate> Name of the radical candidate. It is assumed that
      the files related with this candidate are in a directory with the same
      name.
  -E  <experimental values> guess of the gvalues obtained experimentally in
      python list format, f.e. '[2.0062, 2.0055, 2.0022]'
  -e  <experiment1.dat,exp2.dat,...> experimental field vs absorption files in
      which this candidate can play a role.
  -o  <output.dat='spectrum_wo_hyFiCorr.dat'> dat output file where you want to
      save the field vs spectrum.

  -v  verbose.
  -h  prints this message.

This code should produce:

 - A render of each model called opt.png
 - The computed spectrum obtained by easyspin in a file called spectrum_wo_hyFiCorr.dat
 - A plot of the gvalues of all the candidates called gvalues.png
 - A table of the gvalues in gvalues_table.md
 - A file called vmd_image.md in each one of the directories of the candidates.

This code should be executed in the folder containing the candidate, not in the
candidate directly.
"
exit 0
}

create_vmd_image_md() {
  name=$1;
  subline=$(printf "%0.s=" $(seq 1 ${#name}) );
  cat <<EOF > ./vmd_image.md
$name
$subline

<div align="center">
  <img src="./opt.png"  width="500">
</div>

EOF
}

# ----- set up starts ---------------------------------------------------------
# General variables
output="spectrum_wo_hyFiCorr.dat"
experiment=''
verbose='false'
exper_values=''
mol_cand=''
while getopts 'c:E:e:fO:o:vh' flag;
do
  case "${flag}" in
    c) mol_cand=${OPTARG} ;;
    E) exper_values=${OPTARG} ;;
    e) experiment=${OPTARG} ;;
    o) output=${OPTARG} ;;

    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

if [ ! -d $mol_cand ]
then
  fail "You have to provide a the name of a molecule that is a candidate. Use
        the flag -c for this proporsal"
fi

if [ ! -f $experiment ]
then
  fail "You have to give the dat file of the field of the experiment using the
    flag -e. Check 'myutils extract_EPRspec -h' for details."
fi

source "$(myutils basics -path)" CandInfo $verbose

# starting information
verbose "JOB information"
echo " * Date:"
date
echo " * Command:"
echo "$0" "$@"

# ---- BODY --------------- ----------------------------------------------------
reference=$(myutils gval_workflow -path)
reference=${reference%/*}

# Creates gvalues_table.md and gvalues.png
verbose "create gvalues_table.md and gvalues.png"
myutils extract_system_info "$mol_cand" "$exper_values"
cd $mol_cand

# Radicals
for candidate in *-*/
do
  cd $candidate

  verbose "Render image of the candidate $candidate"
  cp $reference/create_mol_png.tcl .
  vmd -e create_mol_png.tcl -args model_opt.xyz opt.png
  rm create_mol_png.tcl
  
  if [ ! -f $output ]
  then
    fail "$output does not exist in $(pwd). Execute 'myutils extract_EPRspec'
      first."
  fi
  
  # Create vmd_image.md
  name=$mol_cand/$candidate
  create_vmd_image_md $name

  mapfile -t all_exper < <(echo -e "${experiment//,/\\n}")
  # add fit to experimental spectrum
  for exper_spect in ${all_exper[@]}
  do
    name_w_ext=${exper_spect##*/}
    name=${name_w_ext%.dat}.png
    myutils spect_w_experiment $exper_spect $output $name
    cat <<EOF >> ./vmd_image.md

<div align="center">
  <img src="./$name"  width="500">
</div>
EOF
  done
  cd ../
done

# Visualization of molecules where frequencies where computed but not radicals

cd ../ # goes completely out, to where the molecules are.

mapfile -t frequencies < <(find . -name 'model_freq.out' | sort )

ori=$(pwd)
for freq in ${frequencies[@]};
do
  if [[ "$freq" != *"-"* ]];
  then
    name=${freq%/*}

    cd $name
    verbose "Render image of the non-rad molecule $name"
    cp $reference/create_mol_png.tcl .
    vmd -e create_mol_png.tcl -args model_opt.xyz opt.png
    rm create_mol_png.tcl

    # Create vmd_image.md
    create_vmd_image_md $name
    cd $ori
  fi
done

finish "finished"