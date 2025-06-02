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
  -e  <file.dat> experimental field vs absorption file.
  -f  Use this flag to take into account hyperfine corrections. This uses a lot
      of RAM memory. Be sure that you have enough memory or that you filtered.
  -O  <file.out='epr_info.out'> orca output file with computed EPR quantities.
  -o  <output.dat='spectrum_wo_hyFiCorr.dat'> dat output file where you want to
      save the field vs spectrum.
  -m  <float=179.813> experimental value of the microwave frequency. The
      default value corresponds to G-band experiments.
  -n  <int=401> number of data points used to predict the absorption spectrum

  -v  verbose.
  -h  prints this message.

This code should produce:

 - A render of each model called opt.png
 - The computed spectrum obtained by easyspin in a file called spectrum_wo_hyFiCorr.dat
 - A plot of the gvalues of all the candidates called gvalues.png
 - A table of the gvalues in gvalues_table.md
 - A file called vmd_image.md in each one of the directories of the candidates.

This code should be executed in the folder containing the candidate.
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
orca_output="epr_info.out"
MicroWaveExper=179.813
ndpoints=401
output="spectrum_wo_hyFiCorr.dat"
experiment=''
verbose='false'
exper_values=''
hyperfine=''
mol_cand=''
while getopts 'c:E:e:fO:o:m:n:vh' flag;
do
  case "${flag}" in
    c) mol_cand=${OPTARG} ;;
    E) exper_values=${OPTARG} ;;
    e) experiment=${OPTARG} ;;
    f) hyperfine='-f' ;;
    O) orca_output=${OPTARG} ;;
    o) output=${OPTARG} ;;
    m) MicroWaveExper=${OPTARG} ;;
    n) ndpoints=${OPTARG} ;;

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
  vmd -e create_mol_png.tcl -args opt.xyz opt.png
  rm create_mol_png.tcl
  
  verbose "Compute the spectrum with easyspin"
  $(myutils extract_EPRspec -path) -e "$experiment" $hyperfine -O $orca_output \
                                   -o $output -m $MicroWaveExper -n $ndpoints \
                                   -v || fail \
                                   "error extracting spectrum of $candidate"
  
  # Create vmd_image.md
  name=$mol_cand/$candidate
  create_vmd_image_md $name
  for exper_spect in ${experiment%/*}/*.dat
  do
    # add fit to experimental spectrum
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

mapfile -t frequencies < <(find . -name 'freq.out' | sort )

ori=$(pwd)
for freq in ${frequencies[@]};
do
  if [[ "$freq" != *"-"* ]];
  then
    name=${freq%/*}

    cd $name
    verbose "Render image of the non-rad molecule $name"
    cp $reference/create_mol_png.tcl .
    vmd -e create_mol_png.tcl -args opt.xyz opt.png
    rm create_mol_png.tcl
  
    # Create vmd_image.md
    create_vmd_image_md $name
    cd $ori
  fi
done
