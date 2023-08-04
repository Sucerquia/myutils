#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
source $(myutils basics -path) TestChecker

print_help() {
echo "
Check that all scripts are already tested. So far, only python and bash scripts
are included. The output reports which functions/methods and are not tested or
if the test file doesn't even exist. This is done assuming that the test of an
<x_script.ext> is called <test_x_script.py>; same for methods and classes.

    -n    <name> Package you want to check. there must be a '<name> path'
          command that returns the location of the pkg files.
    -t    <test_directory> Default: test

    -h    prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
test_directory='tests'

while getopts 'd:f:n:t:h' flag;
do
    case "${flag}" in
      d) ign_dir=${OPTARG} ;;
      f) ign_fil=${OPTARG} ;;
      n) pkg_name=${OPTARG} ;;
      t) test_directory=${OPTARG} ;;

      h) print_help
    esac
done

mod_path=$( $pkg_name path )   # path to the files to be documented


if [ ${#pkg_name} -eq 0 ];
then
  fail "You have to provide the name of the package you want to check"
fi

if [ ! -d $mod_path/$test_directory ];
then
  fail "$mod_path/$test_directory does not exist"
fi

# directories to be ignore during the check.
ignore_dirs=$(echo $ign_dir | tr ',' '\n')
# files to be ignore during the check.
ignore_files=$(echo $ign_fil | tr ',' '\n')


# ==== Directories ============================================================
# search directories in the package
cd $mod_path
for ign_dir in ${ignore_dirs[@]}
do
    bool_ign="$bool_ign  -path '*$ign_dir*' -o"
done
# the names of subdirectories in the tests respect to the package directory
# will have the same structure
pck_dirs=$( eval "find . -type d -not \(" "${bool_ign::-2}" \
                                   "-prune \)" )

for ign_fil in ${ignore_files[@]}
do
    bool_ign="$bool_ign -name '*$ign_fil*' -o"
done
pck_fils=$( eval "find . -type f -not \(" "${bool_ign::-2}" \
                                   "-prune \)" )

with_test=( )
without_test=( )
cd tests
for fil in $pck_fils
do
    ext=$(echo $fil | cut -d '.' -f3 )
    if [ "$ext" == 'py' ];
    then
        [ -f ${fil%/*}/test_${fil##*/} ] && \
            with_test+=( $fil ) || \
            without_test+=( $fil )
    elif [ "$ext" == 'sh' ]
    then
        namefile=${fil##*/}
        namewoext=$( echo $namefile | cut -d '.' -f 1 )
        # search function
        exist_fun=$( grep "def test_$namewoext("  ${fil%/*}/test_*.py | wc -l )
        # second option in the next "or" searches test file
        [ $exist_fun -eq 1 ] || [ -f ${fil%/*}/test_$namewoext.py ] && \
            with_test+=( $fil ) || \
            without_test+=( $fil )
    fi
done

adjust "python functions"

for fil in ${with_test[@]}
do
    functions=$( grep "def " ../$fil | sed "s/ //g" | grep "^def" | \
                 cut -c4- | cut -d "(" -f 1 | grep -v "^_" )
    lacked_funcs=( )
    for func in ${functions[@]}
    do
        if [[ ! $func == "_"* ]]
        then
            [[ $( grep "def test_$func"  ${fil%/*}/test_${fil##*/} | wc -l ) -eq 0 ]] && \
                lacked_funcs+=( $func )
        fi
    done
    if [ ${#lacked_funcs} -ne 0 ]
    then
        name=$( echo $fil | sed "s/^.\//$pkg_name\//g" )
        echo "$name"
        for func in ${lacked_funcs[@]}
        do
            echo "  $func"
        done
    fi
done

echo
adjust "python classes"

for fil in ${with_test[@]}
do
    functions=$( grep "class " ../$fil | sed "s/ //g" | grep "^class" | \
                 cut -c6- | cut -d ":" -f 1 | cut -d "(" -f 1 | grep -v "^_" )
    lacked_funcs=( )
    for func in ${functions[@]}
    do
        if [[ ! $func == "_"* ]]
        then
            [[ $( grep "def test_$func"  ${fil%/*}/test_${fil##*/} | wc -l ) -eq 0 ]] && \
                lacked_funcs+=( $func )
        fi
    done
    if [ ${#lacked_funcs} -ne 0 ]
    then
        name=$( echo $fil | sed "s/^.\//$pkg_name\//g" )
        echo "$name"
        for func in ${lacked_funcs[@]}
        do
            echo "  $func"
        done
    fi
done

echo
adjust "NEXT FILES DO NOT HAVE A TEST (.\_/.)"

for fil in ${without_test[@]}
do
    name=$( echo $fil | sed "s/^.\//$pkg_name\//g" )
    echo $name
done
