#!/bin/bash

# checks that it has the basic structure
source $(myutils basics -path) just_check 'true'

for file in $@
do
  if [ $file  == "./delete_trial.sh" ]
  then
    continue
  fi
  inhelp='false'
  incase='false'
  inflags='false'
  defa='false'
  whencalling='false'

  if grep -q "while getopts" $file
  then
    grep -E -B 1 "h+[[:space:]]+prints this message" $file | head -n 1 | grep -qE "v+[[:space:]]+verbose" || inhelp='true'
    grep -q "vh' flag" $file || incase='true'
    grep -q "verbose='true" $file || inflags='true'
    grep -q "^verbose=" $file || defa='true'
    grep "source" $file | grep "myutils basics" | grep -q "\$verbose" || whencalling='true'

    if [[ "$inhelp" == 'true' ]]
    then
      echo $file inhelp
      
      n=$(grep -n -E -B 1 "h+[[:space:]]+prints this message" $file | head -n 1 | cut -d "-" -f 1)
      if [ ${#n} -eq 0 ]
      then
        echo "OCHAS inhelp in $file"
      else
        n=$(( n + 1 ))
        sed -i "${n}s/^/  -v  verbose.\n/" $file
      fi
    fi

    if [[ "$incase" == 'true' ]]
    then
      echo $file incase
      sed -i "s/h' flag/vh' flag/g" $file
    fi

    if [[ "$inflags" == 'true' ]]
    then
      echo $file inflags
      n=$(grep -n "h) print_help" $file | cut -d ":" -f 1)
      if [ ${#n} -eq 0 ]
      then
        echo "OCHAS h) in $file"
      else
        sed -i "${n}s/^/  v)  verbose='true' ;;\n/" $file
      fi
    fi

    if [[ "$defa" == 'true' ]]
    then
      echo $file defa
      n=$(grep -n "while getopts" $file | cut -d ":" -f 1)
      if [ ${#n} -eq 0 ]
      then
        echo "OCHAS defa $file"
      else
        sed -i "${n}s/^/verbose='false'\n/" $file
      fi
    fi

    if [[ "$whencalling" == 'true' ]]
    then
      echo $file whencalling
      if [ ${#n} -eq 0 ]
      then
        echo "OCHAS whenca $file"
      else
        n=$(grep -n "myutils basics" $file | cut -d ":" -f 1)
      fi
      sed -i "${n}s/$/ \$verbose/" $file
    fi
  fi
done

# TODO: Add checking Job information in files like command and so: checl
# TODO: compute_forces for an example

# TODO: change headings in bash scripts (functions, settings, body)

# TODO: add finish all codes that sources myutils basics
