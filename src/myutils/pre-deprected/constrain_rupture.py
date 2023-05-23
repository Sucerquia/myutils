from myutils.miscellaneous import output_terminal
import sys
import numpy as np


# DEPRECTED
def extract_bonds(logfile):
    init = output_terminal(f'grep -ni "initial parameters" {logfile} | ' +
                        'head -n 1 | cut -d ":" -f 1 ')
    end = output_terminal(f'tail -n +$(( {init} + 5 )) {logfile} | ' +
                        'grep -n -i "^ ---" | head -n 1 | cut -d ":" -f 1 ')
    table_dofs = output_terminal(f'tail -n +$(( {init} + 5 )) {logfile} | ' +
                                f'head -n $(( {end} - 1 )) | grep "! R" | ' +
                                "awk '{ if ($3){print $3}}' | sed 's/R//g'")
    dofs = table_dofs.split('\n')[:-1]

    lengths = []
    for dof in dofs:
        lengths.append(eval(dof))

    lengths = np.array(lengths)

    return lengths

# DEPRECTED
def compare_dofs(log1, log2):
    bonds1 = extract_bonds(log1).tolist()
    bonds2 = extract_bonds(log2).tolist()
    with open('file1.dat', 'w') as file1:
        for bond in bonds1:
            file1.write(str(bond) + '\n')
    with open('file2.dat', 'w') as file2:
        for bond in bonds2:
            file2.write(str(bond) + '\n')
    print(len(bonds1), len(bonds2))
    print("--------------")

    extra_dofs = []
    for bond in bonds2:
        print(bond, end=' ')
        #if ((bond) not in list(bonds1)) or (list(bond[::-1]) not in list(bonds1)):
        if (bond not in bonds1) and (bond[::-1] not in bonds1):
            print(bond[::-1], end=' ')
            extra_dofs.append(bond)
            print('does not belong')
        else:
            print('<----belong')
    return extra_dofs

print('result   ', compare_dofs(sys.argv[1], sys.argv[2]))
