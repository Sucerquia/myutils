#!/usr/bin/bash

# TODO: check if I have to add a verbose in the next line
source "$(myutils basics -path)" test-basics;
"$@"
