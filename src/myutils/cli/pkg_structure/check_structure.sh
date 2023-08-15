#!/bin/bash

print_help() {
echo "
Check the structure of myutils. All checkers run by default.

    -p    pep8 convention in all python scripts.
    -s    ShellCheck in all bash scripts.
    -t    check tests.

    -h    prints this message.
"
exit 0
}

# ----- definition of functions finishes --------------------------------------

all="true"
pep8="false"
shellcheck="false"
tests="false"

while getopts 'psth' flag; do
    case "${flag}" in
      p) pep8='true' ; all="false" ;;
      s) shellcheck='true' ; all="false" ;;
      t) tests='true' ; all="false" ;;

      h) print_help ;;
      *) echo "for usage check: myutils <function> -h" >&2 ; exit 1 ;;
    esac
done

if $all
then
    pep8="true"
    shellcheck="true"
    tests="true"
fi

if $tests
then
    echo ; echo
    # Ignore directories and files for tests checker
    ign_dirs='pycache,cli,examples,doc_scripts,pre-deprected,tutorials,tests'
    ign_fils='__init__.'
    myutils check_tests -n myutils -d $ign_dirs -f $ign_fils
fi

if $pep8
then
    source "$(myutils basics -path)" PEP8
    original_cs=$(pwd)
    cd "$(myutils path)" || fail "myutils path does not exist"
    echo ; echo
    adjust "Starts"
    pycodestyle . --exclude=pre-deprected --ignore W605
    cd "$original_cs" || fail "returning to former directory"
    finish
fi

if $shellcheck
then
    echo ; echo
    myutils bash_style
fi