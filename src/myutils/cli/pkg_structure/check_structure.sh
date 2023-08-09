#!/bin/bash

source $(myutils basics -path) STRUCTURE_CHECKER

# Ignore directories for tests checker
ign_dirs='pycache,cli,examples,doc_scripts,pre-deprected,tutorials,tests'
ign_fils='__init__.'


myutils check_tests -n myutils -d $ign_dirs -f $ign_fils
