#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 molecule.pdb"
    exit 1
fi

PDB=$1
TMPFILE=$(mktemp)

cat > $TMPFILE << EOF
mol new $PDB

set outfile [open bonds_from_vmd.dat w]
set sel [atomselect top "all"]
foreach {i} [\$sel getbonds] {
    puts \$outfile "\$i"
    puts  "\$i"
}

close \$outfile
quit
EOF

vmd -dispdev text -e $TMPFILE
rm $TMPFILE