#!/bin/bash

# ----- definition of functions -----------------------------------------------
print_help() {
echo "
Takes all the g09 log files given as arguments. If the log file does not report
a proper termination, a new job is resubmitted (creating a backup first) using
'myutils opt_and_forces'.

  -v  verbose
  -h  prints this message.
"
exit 0
}

# ----- set up starts ---------------------------------------------------------
# General variables
verbose='false'
while getopts 'vh' flag;
do
  case "${flag}" in
    v) verbose='true' ;;
    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

# ---- BODY -------------------------------------------------------------------

source "$(myutils basics -path)" recover 'true'
original_path=$(pwd)

files=$@
if [ ${#files} -eq 0 ]
then
    files=$(find . -name *-opt.log | sort)
fi

for logfile in ${files[@]}
do
  if [ ! $(grep "Normal termination" $logfile) ]
  then
    # if there is a path:
    if echo "$logfile" | grep -q "/"
    then
      path=${logfile%/*}
    else
      path="./"
    fi
    just_name=${logfile##*/}
    echo $path $just_name
    cd $path
    cp $just_name tmp
    frozen=$(tail -n 1 ${just_name%.log}.com)
    create_bck ${just_name%.log}.*
    create_bck ${just_name%-opt.log}.xyz
    mv tmp $just_name

    myutils log2xyz "$just_name" || fail "Extracting xyz from logfile"
    file=${just_name//-opt.log/.xyz}
    echo $file
    mv ${just_name%.log}.xyz $file
    myutils change_distance \
      $file ${file%.xyz}-opt \
      "nofile" 0 0 "scale_distance" \
      || fail "Preparating g09 input"
    comfile=${file%.xyz}-opt.com
    echo $comfile
    sed -i '$d' $comfile
    echo $frozen >> $comfile
    sed -i "1a %NProcShared=8" "$comfile"
    sed -i "3a opt(modredun,calcfc)" "$comfile"
    if [[ "$(whoami)" == "hits_"* ]]
    then
      single_part="--partition=single"
    else
      single_part=""
    fi
    sbatch --job-name="${file:0:6}_opt" $single_part \
            --output="${file:0:6}_opt.o" \
            --error="${file:0:6}_opt.e" \
      $(myutils opt_and_forces -path) -f ${comfile%.com} -c || \
      fail "submitting Job {file:0:6}"
    cd $original_path
  fi
done

finish
