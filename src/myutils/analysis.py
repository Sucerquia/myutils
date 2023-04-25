from MDAnalysis.analysis.dihedrals import Ramachandran
import MDAnalysis as ma
from ase.io import read
import numpy as np
from ase.geometry.analysis import Analysis
from ase.io import read


def ramachandran(sith, pdb_file):
    """
    Returns the ramachandran angles of several stretched configurations.
    each compoe

    Parameters
    ==========
    sith:
        sith object.
    pdb_file:
        pdb file that contains the structural information. It doesn't matter if
        the coordinates don't correspond to something relevant. The importance
        are in the structural description (aminoacids, number of residue...)
    """
    u = ma.Universe(pdb_file)
    r = Ramachandran(u.select_atoms('protein')).run()
    rama_angles = []
    u = ma.Universe(pdb_file)
    for conf in sith._deformed:
        u.load_new(conf.atoms.positions)
        r = Ramachandran(u.select_atoms('protein')).run()
        rama_angles.extend(r.results.angles.copy())
    rama_angles = np.array(rama_angles)

    return rama_angles


# add2executable
def distance(file, index1, index2):
    """
    This code checks the distances between two atoms in the last configuration
    of a trajectory file (e.g. *.log file from gaussian)

    Parameters
    ==========

    file: str
        name of the file that contains the trajectory.
    index1: int
        index of the first atom to compute distances
    arg2: int
        index of the second atom to compute distances

    Return
    ======
    (float)  Distance between atoms corresponding with atom with index1 and
    index2.

    Execute from terminal using:
    myutils distance arg1 arg2 arg3

    E.g.
    myutils distance optimization.log 1 20
    """
    index1 = int(index1)
    index2 = int(index2)
    atoms = read(file)
    d = atoms.get_distance(index1, index2)

    return d


def indexes_per_aminoacid(pdb_file):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()
    atoms = [line for line in lines if 'ATOM' in line]
    aminoacids = np.genfromtxt(atoms, usecols=5, dtype=int)
    indexes = np.genfromtxt(atoms, usecols=1, dtype=int)
    atoms_per_aminoacids = {}
    for i in range(1, max(aminoacids) + 1):
        atoms_per_aminoacids[i] = []
    for i in range(len(indexes)):
        atoms_per_aminoacids[aminoacids[i]].append(indexes[i])
    return atoms_per_aminoacids


def all_hydrogen_atoms(mol):
    """"
    find the indexes of all the hydrogen atoms in the peptide from an ASE.Atoms
    object.

    Parameters
    ==========
    mol: string
        ASE.Atoms object to extract the hydrogen indexes.

    Return
    ======
    list of indexes.
    """
    if type(mol) is str:
        mol = read(mol)
    indexes = np.where(mol.get_atomic_numbers() == 1)[0]

    return indexes + 1


def dof_classificator(dofs_indexes, atoms_per_aminoacids):
    list_aminos = {}
    for i in range(1, max(atoms_per_aminoacids.keys()) + 1):
        list_aminos[i] = np.array([], dtype=int)
    for i in range(len(dofs_indexes)):
        for j in atoms_per_aminoacids.keys():
            if np.isin(dofs_indexes[i], atoms_per_aminoacids[j]).all():
                list_aminos[j] = np.append(list_aminos[j], i)
                break
    return list_aminos


def dof_classificator2(dofs_indexes, atoms_per_aminoacids):
    """
    Return all degrees of freedom that includes at least one atom of a residue
    """
    list_aminos = {}
    for i in range(1, max(atoms_per_aminoacids.keys()) + 1):
        list_aminos[i] = np.array([], dtype=int)
    for i in range(len(dofs_indexes)):
        for j in atoms_per_aminoacids.keys():
            if np.isin(dofs_indexes[i], atoms_per_aminoacids[j]).any():
                list_aminos[j] = np.append(list_aminos[j], i)
                break
    return list_aminos


def classical_energies(file):
    "Return the classical energies in Hartrees"
    potential_energy = np.loadtxt(file, usecols=4)
    potential_energy = (potential_energy)*1/2600  # 1Ha=2600kJ/mol
    DOFs_energy = np.loadtxt(file, usecols=[1, 2, 3])*1/2600  # 1Ha=2600kJ/mol

    appr_eDOF = np.sum(DOFs_energy, axis=1)
    appr_eDOF = (appr_eDOF)
    return potential_energy, appr_eDOF


def length_energy(sith, aminos_info, atoms_types):
    """
    return the value of the DOF and the energies as the molecule is stretched

    Parameters
    ==========
    sith: sith object
    aminos_info: dic
        amino_info.amino_info of the requiered amino acid.
    atom_types:
        name of the atoms inside the aminoacid that will be studied,
        example ['CA', 'CB']
    """
    defo = sith._deformed[0]
    try:
        i_ric = defo.dimIndices.index((aminos_info[atoms_types[0]],
                                       aminos_info[atoms_types[1]]))
    except ValueError:
        i_ric = defo.dimIndices.index((aminos_info[atoms_types[1]],
                                       aminos_info[atoms_types[0]]))
    energies = sith.energies.T[i_ric]
    values_dof = []
    for defo in sith._deformed:
        values_dof.append(defo.ric[i_ric])
    values_dof = np.array(values_dof)
    return [values_dof, energies]


def le_same_aminoacids(sith, aminos_info, atoms_types, kind_aminoacid):
    indexes = []
    for j, amino_name in enumerate(aminos_info.amino_name.values()):
        if amino_name in kind_aminoacid:
            indexes.append(j+1)
    all_le = []
    for index in indexes:
        values = length_energy(sith, aminos_info.amino_info[index],
                               atoms_types)
        all_le.append(values)
    return all_le


# ------------------ remove ---------------------------------------------------
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


def extract_bonds(readable_file):
    atoms = read(readable_file)
    ana = Analysis(atoms)
    bonds_zm = ana.unique_bonds

    bonds = []
    for i, pairs in enumerate(bonds_zm[0]):
        if len(pairs) != 0:
            for pair in pairs:
                bonds.append([i+1, pair+1])
    return bonds


# add2executable
def diff_bonds(file1, file2, file='frozen_dofs.dat'):
    """
    This function returns the bonds that are in one file (first argument) but
    not in the other (second argument).
    """

    bonds1 = extract_bonds(file1)
    bonds2 = extract_bonds(file2)

    different_bonds = []
    for bond in bonds1:
        if bond not in bonds2:
            different_bonds.append(bond)

    with open(file, 'a') as fdofs:
        for bond in different_bonds:
            for index in bond:
                fdofs.write(str(index) + ' ')
            fdofs.write('F\n')
    return different_bonds
