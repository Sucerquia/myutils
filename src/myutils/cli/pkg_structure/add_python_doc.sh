#!/bin/bash


print_help() {
echo "
Code that explores the files in the package and automatically creates the
documentation of all classes and functions that finds in it.

   -d   <dir1,dir2...> directories to be ignored.
        Default: 'pycache,tests,tutorials,pre-deprected'
   -f   <fil1,fil2...> files to be ignored. Default: '__init__'
   -p   <absolute_path> path directory to be checked (no relative path).
        Default: \"\$myutils path\"
   -n   <name> pkg name. Default: myutils 

   -v   verbose.
   -h   prints this message.
"
exit 0
}

insert_doc() {
  doc_num_start=$1
  fil=$2
  func=$3
  # Remove prev documentation first
  if awk -v numline=$doc_num_start 'NR==numline' \
          $fil | grep -q "\"\"\""
  then
    doc_num_end=$(tail -n +$(( doc_num_start + 1 )) $fil | \
                  grep -n "\"\"\"" | head -n 1 | cut -d ":" -f 1)
    doc_num_end=$(( doc_num_end + doc_num_start ))
    sed -i "${doc_num_start},${doc_num_end}d" $fil
  fi

  # insert new documentation
  sed -i "$(( doc_num_start - 1 ))r final_$func_doc.txt" $fil

  # delete documentation file
  rm final_$func_doc.txt || fail "not final documentation found final_$func_doc.txt"
}
# ----- definition of functions finishes --------------------------------------

# ==== General variables ======================================================
mod_path=""   # path to the dir with the files to be documented
# directories to be ignored during documentation.
raw_ign_dirs='pycache,tests,ipynb_checkpoints,tutorials,pre-deprected'
# files to be ignored during the documentation.
raw_ign_fils='__init__.'
pkg_name="myutils"

# ==== Costumer set up ========================================================
directory="$(myutils path)"
verbose='false'
while getopts 'd:f:m:n:p:vh' flag;
do
    case "${flag}" in
      d) raw_ign_dirs=${OPTARG} ;;
      f) raw_ign_fils=${OPTARG} ;;
      n) pkg_name=${OPTARG};;
      p) mod_path=${OPTARG};;

      v) verbose='true' ;;
      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done
source "$(myutils basics -path)" AddPythonDoc $verbose

if [ "$mod_path" == "" ]
then
  mod_path=$(myutils path)
fi

mapfile -t ignore_dirs < <(echo "$raw_ign_dirs" | tr ',' '\n')
# files to be ignore during the check.
mapfile -t ignore_files < <(echo "$raw_ign_fils,$raw_ign_dirs" | tr ',' '\n')

# ==== Directories ============================================================
cd "$mod_path" || fail "$mod_path not found"

# directories to ignore
for ign_dir in "${ignore_dirs[@]}"
do
  bool_ign="$bool_ign -path '*$ign_dir*' -o"
done
# files to ignore
for ign_fil in "${ignore_files[@]}"
do
    bool_ign="$bool_ign -name '*$ign_fil*' -o"
done

# files
mapfile -t pck_fils < <(eval "find . -type f -not \(" "${bool_ign::-2}" \
                        "-prune \)" )

for fil in "${pck_fils[@]}"
do
  verbose $fil
  ext=$(echo "$fil" | cut -d '.' -f3 )
  # Python files
  if [ "$ext" == 'py' ]
  then
    module=$(echo "myutils"${fil//\.\//\.} | sed "s/\//\./g" | sed "s/\.py//g")
    # Functions
    mapfile -t functions < <(grep "^def " "$fil" | awk '{print $2}' | \
                             cut -d "(" -f 1)
    if [ ${#functions} -eq 0 ]; then verbose "functions"; fi

    for func in ${functions[@]}
    do
      echo $func
      # The output of the next function is stored in final_<func>_doc.txt
      myutils python_doc_fixer -f "$func" -m "$module" || \
        fail "creating new documentation"

      # search n lines of the beginning of the function, the end of the
      # heading of the function and the beginning of the documentation
      n_func=$( grep -n "def $func" $fil | cut -d ":" -f 1 )
      rel_n_func_end=$( tail -n +$n_func $fil | grep -n ")" | head -n 1 | \
                        cut -d ':' -f 1)
      doc_num_start=$(( n_func + rel_n_func_end ))

      insert_doc $doc_num_start $fil $func
    done

    # Classes
    mapfile -t classes < <(grep "^class " "$fil" | awk '{print $2}' | \
                             cut -d "(" -f 1)
    if [ ${#classes} -eq 0 ]; then verbose "Classes"; fi
    for class in ${classes[@]}
    do
      # Documentation of the class definition
      class=${class%:}
      echo $class
      myutils python_doc_fixer -m $module -c $class || \
        fail "creating new documentation of $class"

      # search n lines of the beginning of the function, the end of the
      # heading of the function and the beginning of the documentation
      n_class=$( grep -n "^class $class" $fil | cut -d ":" -f 1 )
      rel_n_class_end=$( tail -n +$n_class $fil | grep -n ":" | head -n 1 | \
                        cut -d ':' -f 1)
      doc_num_start=$(( n_class + rel_n_class_end ))

      insert_doc $doc_num_start $fil

      # TODO: add description of the attributes.

      # ==== description of modules
      # subfile with the class only
      n_end_class=$(tail -n +$n_class $fil | grep -nEv "^ |^$"  | \
                    grep -v "^$" | tail -n +2 | head -n 1 | cut -d ":" -f 1 )

      if [ ${#n_end_class} -eq 0 ]
      then
        tail -n +$n_class $fil > ${class}_complete.dat
      else
        myutils find_blocks -f $fil -s $n_class -e $n_end_class \
                            -i -o ${class}_complete
      fi
      mapfile -t methods < <(myutils methods_in_class $module $class | \
                             head -n -1)
      echo "$class methods"
      for method in ${methods[@]}
      do
        echo $method
        myutils python_doc_fixer -m $module -c $class -f $method || \
          fail "creating new documentation of $class"
        rel_n_meth=$(grep -n "def $method" ${class}_complete.dat)
        rel_n_meth_end=$( tail -n +$rel_n_meth  | grep -n ")" | head -n 1 | \
                          cut -d ':' -f 1)
        doc_num_start=$(( n_class + rel_n_meth + rel_n_meth_end ))
        insert_doc $doc_num_start $fil $method
      done
      rm ${class}_complete.dat
    done
  fi
done

finish
