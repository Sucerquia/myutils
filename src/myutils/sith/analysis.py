import numpy as np


def dof_classificator_all(dofs_indexes, atoms_per_aminoacids):
    """
    Separates all degrees of freedom defined by atoms (all) of a residue.

    Parameters
    ==========
    dof_indexes: list of duples
        sith._deformed[n].dimIndices containing definition of the degrees of
        freedom in term of the atomic indexes.
    atoms_per_aminoacids: dict
        Atoms in each residue. The keys are the number of the residues, values
        should be the indexes of the atoms belonging to the residue of the key.

    Return
    ======
    (dict) [keys: Residues (int), values: (list) [#DOFsPerResidue (int)]]
    Indexes of the degrees of freedom containing all atoms of each residue.

    Note
    ====
    atoms_per_aminoacids can be obtained from
    myutils.peptides.PepSetter.atom_indexes
    """
    list_aminos = {}
    for i in range(1, max(atoms_per_aminoacids.keys()) + 1):
        list_aminos[i] = np.array([], dtype=int)
    for i in range(len(dofs_indexes)):
        for j in atoms_per_aminoacids.keys():
            if np.isin(dofs_indexes[i], atoms_per_aminoacids[j]).all():
                list_aminos[j] = np.append(list_aminos[j], i)
                break
    return list_aminos


def dof_classificator_one(dofs_indexes, atoms_per_aminoacids):
    """
    Separates all degrees of freedom defined by atoms (at least one) of a
    residue.

    Parameters
    ==========
    dof_indexes: list of duples
        sith._deformed[n].dimIndices containing definition of the degrees of
        freedom in term of the atomic indexes.
    atoms_per_aminoacids: dict
        Atoms in each residue. The keys are the number of the residues, values
        should be the indexes of the atoms belonging to the residue of the key.

    Return
    ======
    (dict) [keys: Residues (str), values: (list) [#DOFsPerResidue (int)]]
    Indexes of the degrees of freedom containing at least one atom of each
    residue.

    Note
    ====
    atoms_per_aminoacids can be obtained from
    myutils.peptides.PepSetter.atom_indexes
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


def length_energy(sith, aminos_info, atom_types):
    """
    Return distances between two atom types in one amino acid and the
    energy associated with this DOF as the molecule is stretched.

    Parameters
    ==========
    sith: sith object
        sith object containing the distribution of energies. That implies to
        have the class variable 'energies' with the energies per deformed
        configuration and and per DOF.
    aminos_info: dic
        name of the atoms of one amino acid associated with the index.
    atom_types: str
        name of the atoms inside the aminoacid that will be studied,
        example ['CA', 'CB'].

    Return
    ======
    (list) [2 x #Def (float)] values of the DOF and energies associate with
    those DOFs per deformed configuration in the selected amino.

    Note
    ====
    aminos_info can be obtained from
    myutils.peptides.PepSetter.amino_info[n] where n is the selected amino
    acid.
    """
    defo = sith._deformed[0]
    try:
        i_ric = defo.dimIndices.index((aminos_info[atom_types[0]],
                                       aminos_info[atom_types[1]]))
    except ValueError:
        i_ric = defo.dimIndices.index((aminos_info[atom_types[1]],
                                       aminos_info[atom_types[0]]))
    energies = sith.energies.T[i_ric]
    values_dof = []
    for defo in sith._deformed:
        values_dof.append(defo.ric[i_ric])
    values_dof = np.array(values_dof)
    return [values_dof, energies]


def le_same_aminoacids(sith, peptides_info, atom_types, kind_amino):
    """
    Return distances between two atom types in the same type of amino acid and
    the energy associated with these DOF as the molecule is stretched.

    Parameters
    ==========
    sith: sith object
        sith object containing the distribution of energies. That implies to
        have the class variable 'energies' with the energies per deformed
        configuration and and per DOF.
    peptides_info:
        object with the info of the peptide.
    atom_types: str
        name of the atoms inside the aminoacid that will be studied,
        example ['CA', 'CB'].
    kind_amino:
        name of the amino acid.

    Return
    ======
    (list) [#kindAmino x [2 x #Def (float)]] values of the DOF and energies
    associate with those DOFs per deformed configuration in the selected amino.

    Note
    ====
    peptides_info can be obtained from
    myutils.peptides.PepSetter
    """
    indexes = []
    for j, amino_name in enumerate(peptides_info.amino_name.values()):
        if amino_name in kind_amino:
            indexes.append(j+1)
    all_le = []
    for index in indexes:
        values = length_energy(sith, peptides_info.amino_info[index],
                               atom_types)
        all_le.append(values)
    return all_le
