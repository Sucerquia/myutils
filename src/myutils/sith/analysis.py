import numpy as np
import matplotlib.pyplot as plt
from myutils.plotters import StandardPlotter
from SITH.SITH import SITH
from myutils.miscellaneous import output_terminal
from myutils.peptides import PepSetter
from os.path import isdir


def dof_classificator_all(dofs_indexes, atoms_per_aminoacids):
    """
    Separates all degrees of freedom defined by atoms (all) of a residue.

    Parameters
    ==========
    dof_indexes: list of duples
        sith.structures[n].dim_indices containing definition of the degrees of
        freedom in term of the atomic indexes.
    atoms_per_aminoacids: dict
        Atoms in each residue. The keys are the number of the residues, values
        should be the indexes of the atoms belonging to the residue of the key.

    Return
    ======
    (dict) [keys: Residues (int), values: (array) [#DOFsPerResidue (int)]]
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
        sith.structures[n].dim_indices containing definition of the degrees of
        freedom in term of the atomic indexes.
    atoms_per_aminoacids: dict
        Atoms in each residue. The keys are the number of the residues, values
        should be the indexes of the atoms belonging to the residue of the key.

    Return
    ======
    (dict) [keys: Residues (str), values: (array) [#DOFsPerResidue (int)]]
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
    defo = sith.structures[0]
    try:
        i_ric = defo.dim_indices.index((aminos_info[atom_types[0]],
                                       aminos_info[atom_types[1]]))
    except ValueError:
        i_ric = defo.dim_indices.index((aminos_info[atom_types[1]],
                                       aminos_info[atom_types[0]]))
    energies = sith.energies.T[i_ric]
    values_dof = []
    for defo in sith.structures:
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
    kind_amino: (list) [(str)]
        name of the amino acids.

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
            indexes.append(j + 1)
    all_le = []
    for index in indexes:
        values = length_energy(sith, peptides_info.amino_info[index],
                               atom_types)
        all_le.append(values)
    return all_le


# TODO: deprecated, up to date with last version of sith.
class SithAnalysis:
    def __init__(self, sith, pepinfo):
        self.sith = sith
        self.pep_info = pepinfo

    def le_dof_amino(self, a_names, aminos):
        if isinstance(aminos, int):
            # if all atoms belog to the same aminoacid
            aminos = [aminos for _ in a_names]
        else:
            assert len(aminos) == len(a_names)

        indexes = []
        for amino, atom in zip(aminos, a_names):
            indexes.append(self.pep_info.amino_info[amino][atom])

        dof = np.zeros(4, dtype=int)
        dof[-len(indexes):] = indexes
        dof_i = self.index_dof(dof)

        dof_value = self.sith.all_dofs[:, dof_i]
        dof_e = self.sith.dofs_energies[:, dof_i]
        dof_e -= dof_e[0]

        return dof_value - dof_value[0], dof_e

    def index_dof(self, target: tuple):
        """
        Search the index of a specific dof.

        Parameter
        =========
        target: np.ndarray
            degree of freedom.

        Return
        ======
        (int) index
        """
        for i, dof in enumerate(self.sith.structures[0].dim_indices):
            dof_wo_0 = dof[np.nonzero(dof)[0]]
            target_wo_0 = target[np.nonzero(target)[0]]

            if (dof_wo_0 == target_wo_0).all() or \
               (dof_wo_0 == target_wo_0[::-1]).all():
                return i
        raise ValueError("Non-found dof.")


class DataSetAnalysis:
    def __init__(self, data_dir='./', exclude_prolines=True, kind_analysis='sith'):
        jedi = False
        sith = False
        if kind_analysis == 'sith':
            sith = True
        elif kind_analysis == 'jedi':
            jedi= True
        elif kind_analysis == 'both':
            jedi = True
            sith = True
        else:
            raise ValueError("Non-recognized kind of analysis.")

        peptides = output_terminal('ls ' + data_dir, print_output=False).split('\n')
        self.pep_infos = []
        self.outcomes = {'jedi': [], 'sith': []}
        self.analysis = {'jedi': [], 'sith': []}

        for pep in peptides:
            print(pep + ' ', end='')
            if ((exclude_prolines) and ('P' in pep)) or (not isdir(f'{data_dir}/{pep}/forces')):
                continue

            self.pep_infos.append(PepSetter(f'{data_dir}/{pep}/{pep}-stretched00.pdb'))
            
            if jedi:
                path = f'{data_dir}/{pep}/'
                sith = SITH(path)
                sith.killer(killElements='H')
                sith.rem_first_last(from_last_minimum=True)
                self.outcomes['jedi'].append(sith)
                self.analysis['jedi'].append(SithAnalysis(self.outcomes['jedi'][-1],
                                                          self.pep_infos[-1]))
            if sith:
                try:
                    path = f'{data_dir}/{pep}/forces'
                    sith = SITH(inputfiles=path)
                    sith.rem_first_last(from_last_minimum=True)
                    sith.sith_analysis()
                    self.outcomes['sith'].append(sith)
                    self.analysis['sith'].append(SithAnalysis(self.outcomes['sith'][-1],
                                                              self.pep_infos[-1]))
                except:
                    print(pep)
                    continue

    def plot_le(self, a_names, aminos=3, ax: plt.Axes = None, sp=None,
                kind_analysis='sith', lw=1, ms=1):
        if sp is None:
            sp = StandardPlotter()
        if ax is None:
            ax = sp.ax[0]
        sp.axis_setter(ax=ax,
                        xlabel=f'Distance({", ".join(a_names)})[A]',
                        ylabel='Energy[Ha]')
        for an in self.analysis[kind_analysis]:
            l, e = an.le_dof_amino(a_names, aminos)
            sp.plot_data(l, e, ax=ax, lw=lw, markersize=ms)
        return ax


def set_hes_from_ref(geo_ref, sith_tar, structure):
    """
    Set the hessian in in a target sith taken from a geometry of reference.
    
    Parameters
    ==========
    geo_ref: SITH.Utilities.Geometry
        sith object that contains the atribute you want to redefine.
    sith_tar: SITH.SITH
        sith object that will change its property.
    structure: int
        index of the deformed structure to set the hessian.

    Returns
    =======
    (SITH.SITH) returns the sith_tar with the hessian in the defined structure.
    
    Note: All the SITH.SITH.structures are Geometry objects with all the information of the structure.
    """
    for dof in sith_tar.dim_indices:
        test = dof[dof != 0]
        check2 = np.concatenate((test[::-1], np.zeros(4-len(test), dtype=int)))
        if not (np.all(geo_ref.dim_indices == dof, axis=1).any() \
            or np.all(geo_ref.dim_indices == check2, axis=1).any()):
            raise('this dof does not exist: ', dof)
        
    order = []
    for dof in sith_tar.dim_indices:
        test = dof[dof != 0]
        check2 = np.concatenate((test[::-1], np.zeros(4-len(test), dtype=int)))
        try:
            index = np.where(np.all(geo_ref.dim_indices == dof, axis=1))[0][0]
        except IndexError:
            index = np.where(np.all(geo_ref.dim_indices == check2, axis=1))[0][0]
        order.append(index)
        
    geo_ref.hessian = geo_ref.hessian[order]
    geo_ref.hessian = geo_ref.hessian[:, order]
    
    sith_tar.structures[structure].hessian = geo_ref.hessian
    
    return sith_tar