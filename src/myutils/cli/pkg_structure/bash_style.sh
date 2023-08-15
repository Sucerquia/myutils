#!/bin/bash

source "$(myutils basics -path)" BASH_CHECKER

bsoriginal_dir=$(pwd)

cd "$(myutils path)" || fail "moving to package directory"

mapfile -t bash_files < <(find . -name "*.sh")

for fil in "${bash_files[@]}"
do
    shellcheck -e SC1090,SC2015 "$fil"
done
verbose

cd "$bsoriginal_dir" || fail "original directory lost"