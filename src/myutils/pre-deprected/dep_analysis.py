from myutils.miscellaneous import output_terminal


def cap_hydrogen_atoms(pdb_file):
    """"
    find the indexes of the hydrogen atoms of the caps parts in the peptide.

    Parameters
    ==========
    pdb_file: string
        name of the pdb file with the molecule of interest.

    Return
    ======
    list of indexes.
    """
    output = output_terminal(f'grep ATOM {pdb_file} | grep -n ATOM | grep -e \
                               ACE -e NME | grep HH')
    output = output.split('\n')[:-1]
    indexes = [int(line.split(':')[0]) for line in output]

    return indexes
