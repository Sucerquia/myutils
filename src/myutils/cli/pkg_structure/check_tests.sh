#!/usr/bin/bash

# ----- definition of functions starts ----------------------------------------
source "$(myutils basics -path)" TestChecker

adjust "Starts"
print_help() {
echo "
Check that all scripts in a package are already tested. So far, only python and
bash scripts are included. The output reports which functions/methods and are
not tested or if the test file doesn't even exist. This is done assuming that
the test of an <x_script.ext> is called <test_x_script.py>; same for methods
and classes.

    -d    <ignore_dirs> directories to be ignored.
    -f    <ignore_files> files to be ignored.
    -p    <dir_path> Path to the package source to be analyzed.
          Default: myutils -path
    -t    <test_directory> subdirectory of the package where tests are located.
          Default: test.

    -h    prints this message.
"
exit 0
}

ignore_elements () {
    # function that removes elements of an array that contains certain keywords
    # n1,n2: number of elements in array 1 and array 2
    # $1: array of elements
    # $2: array of ignore keywords
    mapfile -t ninarr < <(echo "$1" | tr ',' '\n')
    array1=( "${@:2:ninarr[0]}" )
    array2=( "${@:$(( ninarr[0] + 2 ))}" )

    new_array=()
    for element in "${array1[@]}"
    do
        matched=0
        for keyword in "${array2[@]}"
        do
            if [[ "$element" == *"$keyword"* ]]
            then
                matched=1
                break
            fi
        done
        if [[ "$matched" -eq 0 ]]
        then
            new_array+=( "$element" )
        fi
    done
    for element in "${new_array[@]}"; do echo "$element"; done
}


# ----- definition of functions finishes --------------------------------------

# ----- set up starts ---------------------------------------------------------
# General variables
test_directory='tests'
original_cs=$(pwd)
mod_path=$( myutils path )

while getopts 'd:f:p:t:h' flag;
do
    case "${flag}" in
      d) raw_ign_dirs=${OPTARG} ;;
      f) raw_ign_fils=${OPTARG} ;;
      p) mod_path=${OPTARG} ;;
      t) test_directory=${OPTARG} ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

VERBOSE "Tests of $mod_path will be checked"

if [ ! -d "$mod_path/$test_directory" ];
then
  fail "$mod_path/$test_directory does not exist"
fi

# directories to be ignore during the check.
mapfile -t ignore_dirs < <(echo "$raw_ign_dirs" | tr ',' '\n')
# files to be ignore during the check.

mapfile -t ignore_files< <(echo "$raw_ign_fils,$raw_ign_dirs" | tr ',' '\n')

# ==== Directories ============================================================
# search directories in the package
cd "$mod_path" || fail "$mod_path does not exist"


# the names of subdirectories in the tests respect to the package directory
# will have the same structure
mapfile -t all_pkg_dirs < <(find . -type d)

mapfile -t pkg_dirs < <( ignore_elements \
                         "${#all_pkg_dirs[@]},${#ignore_dirs[@]}" \
                         "${all_pkg_dirs[@]}" "${ignore_dirs[@]}" )

mapfile -t all_pkg_files < <(find . -type f)
mapfile -t pkg_files < <( ignore_elements \
                         "${#all_pkg_files[@]},${#ignore_files[@]}" \
                         "${all_pkg_files[@]}" "${ignore_files[@]}" )

with_test=( )
without_test=( )
cd $test_directory

# the structure of package and tests must be the same
adjust "Tests directories structure"
for dir in "${pkg_dirs[@]}"
do
    [ -d "$dir" ] || warning "$dir does not exist in tests"
done

for file in "${pkg_files[@]}"
do
    # Find the existing test files of python scripts
    ext=$(echo "$file" | cut -d '.' -f3 )
    if [ "$ext" == 'py' ]
    then
        [ -f "${file%/*}/test_${file##*/}" ] && \
            with_test+=( "$file" ) || \
            without_test+=( "$file" )
    # Check if which sh files are tested already
    elif [ "$ext" == 'sh' ]
    then
        namefile="${file##*/}"
        namewoext=$( echo "$namefile" | cut -d '.' -f 1 )
        # search function
        # disable SC2126 because it is rearching in several files
        # shellcheck disable=SC2126
        exist_fun=$( grep "def test_$namewoext("  "${file%/*}"/test_*.py  | wc -l )

        # second option in the next "or" searches test file
        [ "$exist_fun" -eq 1 ] || [ -f "${file%/*}/test_$namewoext.py" ] && \
            with_test+=( "$file" ) || \
            without_test+=( "$file" )
    fi
done

adjust "python functions"

for fil in "${with_test[@]}"
do
    mapfile -t functions < <(grep "def " "../$fil" | sed "s/ //g" | \
                             grep "^def" | cut -c4- | cut -d "(" -f 1 | \
                             grep -v "^_" )
    lacked_funcs=( )
    for func in "${functions[@]}"
    do
        if [[ ! "$func" == "_"* ]]
        then
            [[ $( grep -c "def test_$func" "${fil%/*}/test_${fil##*/}" ) \
               -eq 0 ]] && lacked_funcs+=( "$func" )
        fi
    done
    if [ "${#lacked_funcs}" -ne 0 ]
    then
        name=$( echo "$fil" | sed "s/^.\//$pkg_name\//g" )
        for func in "${lacked_funcs[@]}"
        do
            echo "  $func"
        done
    fi
done

adjust "python classes"
for fil in "${with_test[@]}"
do
    mapfile -t classes < <(grep "class " "../$fil" | sed "s/ //g" | \
                           grep "^class" | cut -c6- | cut -d ":" -f 1 | \
                           cut -d "(" -f 1 | grep -v "^_" )
    [ "${#classes}" -eq 0 ]
    lacked_funcs=( )
    for class in "${classes[@]}"
    do
        if [[ ! "$class" == "_"* ]]
        then
            [[ $( grep -c "def test_$class"  "${fil%/*}/test_${fil##*/}" ) -eq 0 ]] && \
                lacked_funcs+=( "$func" )
        fi
    done
    if [ "${#lacked_funcs}" -ne 0 ]
    then
        name=$( echo "$class" | sed "s/^.\//$pkg_name\//g" )
        for func in "${lacked_funcs[@]}"
        do
            echo "  $class"
        done
    fi
done

adjust "NEXT FILES DO NOT HAVE A TEST (.\_/.)"

for fil in "${without_test[@]}"
do
    name=$( echo "$fil" | sed "s/^.\//$pkg_name\//g" )
    echo "$name   "
done

cd $original_cs
