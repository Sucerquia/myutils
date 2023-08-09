#!/bin/bash

source $(myutils basics -path) BASH_CHECKER

original_dir=$(pwd)

cd $(myutils path)

bash_files=$( find . -name "*.sh")
for fil in "${bash_files[@]}"
do
    verbose $fil
    shellcheck -e SC1090 $fil
done
verbose

cd $original_dir