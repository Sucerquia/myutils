#!/bin/bash

source "$(myutils basics -path)" AddPythonDoc

print_help() {
echo "
Code takes the documentation of a function and checks the documentation adding
TODOs in the missing parts. The output is stored in a file called
final_<method>-doc.txt

   -f   <method> Function to be checked
   -m   <module> Module that contains the Function
   -s   <n_spaces> number of leading spaces.

   -h   prints this message.
"
exit 0
}
# ----- definition of functions finishes --------------------------------------

# ==== Costumer set up ========================================================
directory="$(myutils path)"
spaces=0
while getopts 'f:m:s:h' flag;
do
  case "${flag}" in
    f) function=${OPTARG} ;;
    m) module=${OPTARG} ;;
    s) num_spaces=${OPTARG} ;;

    h) print_help ;;
    *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
  esac
done

# ==== Body ===================================================================

# ==== Initial Blocks =========================================================
myutils function_doc $module $function | sed 's/^/#new_line/' > \
  $function-doc.txt
sed -i 's/[[:space:]]*$//g' $function-doc.txt

mapfile -t ns_empty < <(cat $function-doc.txt | grep -n "#new_line$" | \
                        cut -d ":" -f1)

# ==== Existing blocks in old documentation
for (( i=0 ; i < $(( ${#ns_empty[@]} - 1 )) ; i++ ))
do
  myutils find_blocks -f $function-doc.txt \
                      -s ${ns_empty[$i]} \
                      -e ${ns_empty[$(( i + 1 ))]} \
                      -i -o documentation-blocks_$i
done

# ==== Block of Parameters in old documentation
# In case it does not exist, a new block is created
par_block=$(grep -xl "#new_line    Parameters" documentation-blocks_*)
if [ ${#par_block} -eq 0 ]
then
  echo "didn't find the parameters"
  cat << EOF > documentation-blocks_parameters.out
#new_line    Parameters
#new_line    ==========
EOF
  par_block="documentation-blocks_parameters.out"
fi

# ==== Block of Return in old documentation
# In case it does not exist, a new block is created
return_block=$(grep -xl "#new_line    Return" documentation-blocks_*)
if [ ${#return_block} -eq 0 ]
then
  cat << EOF > documentation-blocks_return.out
#new_line    Return
#new_line    ======
EOF
  return_block="documentation-blocks_return.out"
fi

# ==== Block of Definition in old documentation
# In case it does not exist, a new block is created
if [[ "$par_block" == "documentation-blocks_0.out" ]] || \
   [[ "$return_block" == "documentation-blocks_0.out" ]] || \
   [ ! -f "documentation-blocks_0.out" ]
then
  echo "#new_line    # TODO: Add definition" > documentation-blocks_def.out
  definition_block="documentation-blocks_def.out"
else
  definition_block="documentation-blocks_0.out"
fi

# === Check parameters ========================================================
# parameters
mapfile -t parameters < <(myutils args_and_defaults $module $function | \
    grep -v "###" | grep -vx '' )

# insert missed parameters
for par in "${parameters[@]}"
do
  par_name=$(echo $par | cut -d ":" -f 1) # name of the parameter
  par_defa=$(echo $par | cut -d ":" -f 2) # default of the parameter
  n_par=$(grep -n $par_name $par_block | \
          cut -d ":" -f 1) # line number of the parameter

  if [ ${#n_par} -eq 0 ]
  then
    # if the variable is not defined
    echo "#new_line    $par # TODO: check default value" >> $par_block
    echo "#new_line        # TODO: add documentation of this parameter" >> \
      $par_block
  else
    # check if the default exists
    if [ ! ${#par_defa} -eq 0 ]
    then
      awk -v line=$n_par 'NR==line' $par_block | grep -q Default || \
        sed -i "${n_par}s/$/\.$par_defa # TODO: check default value/" $par_block
    fi

    # check if the definition of the parameter exist
    if awk -v line=$(( n_par + 1 )) 'NR==line' $par_block | grep -q ":"
    then
      sed -i "${n_par}a\#new_line        # TODO: add documentation of this parameter" $par_block
    fi
  fi
done

# === Create the final doc block ==============================================
echo "#new_line    \"\"\"" > final_$function-doc.txt
cat $definition_block >> final_$function-doc.txt
rm $definition_block

echo "#new_line" >> final_$function-doc.txt
cat $par_block >> final_$function-doc.txt
rm $par_block

echo "#new_line" >> final_$function-doc.txt
cat $return_block >> final_$function-doc.txt
rm $return_block

# rest of blocks
for other_doc_block in documentation-blocks*.out
do
  echo "#new_line" >> final_$function-doc.txt
  cat $other_doc_block >> final_$function-doc.txt
  rm $other_doc_block
done
echo "#new_line    \"\"\"" >> final_$function-doc.txt

# cleaning
sed -i 's/#new_line//g' final_$function-doc.txt
spaces=$(printf "%${num_spaces}s")
# Add spaces to the beginning of each line and save to a new file
sed -i "s/^/${spaces}/" final_$function-doc.txt
rm $function-doc.txt

finish
