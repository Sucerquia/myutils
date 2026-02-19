#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 molecule.pdb"
    exit 1
fi

PDB=$1
TMPFILE=$(mktemp)

cat > $TMPFILE << EOF
mol new $PDB

set sel [atomselect top "all"]
foreach {i} [\$sel getbonds] {
    puts  "\$i"
}

quit
EOF

vmd -dispdev text -e $TMPFILE
rm $TMPFILE