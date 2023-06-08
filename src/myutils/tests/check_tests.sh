# ==== General variables ======================================================
pkg_name='myutils'
mod_path=$($pkg_name path)   # path to the files to be documented
# directories to be ignore during documentation.
ignore_dirs=( 'pycache' 'tests' 'cli' 'examples' 'doc_scripts' 'pre-deprected')
# files to be ignore during documentation.
ignore_files=( '__init__.' )

source $(myutils basics) TestChecker
# ==== General variables ======================================================

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
    bool_ign="$bool_ign  -name '*$ign_fil*' -o"
done
pck_fils=$( eval "find . -type f -not \(" "${bool_ign::-2}" \
                                   "-prune \)" )

with_test=( )
without_test=( )
cd tests
for fil in $pck_fils
do
    ext=$(echo $fil | cut -d '.' -f3 )
    if [ $ext == 'py' ]
    then
        [ -f ${fil%/*}/test_${fil##*/} ] && \
            with_test+=( $fil ) || \
            without_test+=( $fil )
    elif [ $ext == 'sh' ]
    then
        namefile=${fil##*/}
        namewoext=$( echo $namefile | cut -d '.' -f 1 )
        exist=$( grep "def test_$namewoext("  ${fil%/*}/test_*.py | wc -l )
        [ $exist -eq 1 ] && \
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
