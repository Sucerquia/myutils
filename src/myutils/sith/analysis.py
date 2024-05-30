import numpy as np
import matplotlib.pyplot as plt
from myutils.plotters import StandardPlotter
from SITH.SITH import SITH
from myutils.peptides import PepSetter
from pathlib import Path
from typing import Union, Callable


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

        return dof_value, dof_e

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

            if len(dof_wo_0) != len(target_wo_0):
                continue

            if (dof_wo_0 == target_wo_0).all() or \
               (dof_wo_0 == target_wo_0[::-1]).all():
                return i
        raise ValueError("Non-found dof.")



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


class DataSetAnalysis:
    def __init__(self, inner_steps: Callable,  data_dir: str='./',
                 subdir='',
                 exclude_prolines=True,
                 exclude=None,
                 pdb_pattern='stretched00',
                 struc_pattern='forces*.fchk',
                 **kwargs):
        """
        Creates the objects of to do the analysis of a bunch of peptides in a
        data set.

        Prameters
        =========
        data_dir: str or Path
            directory containing the peptides you want to analyse. HeadUP: this
            directory has to contain a set of directories with the name of all
            the peptides (expected as amino code) and nothing else.
        exclude_prolines: Bool. Default=True
            False to include peptides with prolines.
        pdb_pattern: str. Default='stretched00'
            every peptide directory has to contain one and only one pdb file
            containing the molecular information. This code will look for
            "<data_dir>/<peptide>/*pattern*.pdb".
        struc_pattern: str. Default='forces'
            every peptide directory has to contain a set of files containing
            the structure information requiered by 'sith'. This code will look
            for files named as "<data_dir>/<peptide>/*pattern*". It will be
            assumed that the files are organized alphabetically in the order of
            stretched structure.
        inner_steps: function. Default=None
            function with the operations to apply to every sith object -its
            unique argument- after initialized. Note that it has to have the
            analysis you want to apply. For example,
            ```python
            def inner_steps(sith):
                sith.rem_first_last(from_last_minimum=True)
                sith.sith_analysis()

                return sith
            ```
        **kwargs initializinf sith
        """
        if isinstance(data_dir, (Path, str)):
            path = Path(data_dir)
        assert path.is_dir(), f"{data_dir} does not exist."

        peptides = list(path.glob('*/'))
        assert len(peptides) != 0, f"There are not directories in {data_dir}"

        self.pep_infos = []
        self.names = []
        self.outcomes = []
        self.analysis = []

        print("Log SITH analysis:\n")
        
        prolines = []
        errors = []
        if exclude is None:
            exclude = []

        for pep in peptides:
            if (((exclude_prolines) and ('P' in pep.stem)) or (pep.stem in exclude)):
                prolines.append(pep.stem)
                continue
            assert pep.is_dir(), f"{pep} does not exist."
            
            pdb = list(pep.glob(f'*{pdb_pattern}*.pdb'))
            if len(pdb) == 0:
                raise FileNotFoundError(f"Not '*{pdb_pattern}*.pdb' found in"
                                        f"{str(pep)}")
            elif len(pdb) > 1:
                raise ValueError("There are more than one pdb files with the pattern"
                                 f" {str(pep)}/*{pdb_pattern}*.pdb")

            self.pep_infos.append(PepSetter(pdb[0]))
            
            struc = pep / subdir
            structure_files = list(struc.glob(f'*{struc_pattern}*'))
            structure_files.sort()

            try:
                sith = SITH(inputfiles=structure_files, **kwargs)
                sith = inner_steps(sith)
                sith.name = pep.stem
                self.outcomes.append(sith)
                self.analysis.append(SithAnalysis(self.outcomes[-1],
                                                self.pep_infos[-1]))
            except:
                errors.append(pep.stem)

            print(pep.stem + ' ', end='')
        
        print(f"\n--- A total of {len(self.outcomes)} peptides where added to the "
              "analysis")
    
        if len(prolines) > 0:
            print("--- The next peptides were neglected because they have a "
                  "proline at least:")
            [print(pep + ' ', end='') for pep in prolines]
    
        if len(errors) > 0:
            print("\n--- The next peptides did not woked for some reason. "
                  "Check them individually:")
            [print(pep+ ' ', end='') for pep in prolines]
    
        self.test()
    
    def test(self):
        """check that the first sith has the basic variables. It assumes that
        the rest also has it"""
        assert isinstance(self.outcomes[0].structures_scf_energies, np.ndarray)
        assert isinstance(self.outcomes[0].dims, np.ndarray)
        assert isinstance(self.outcomes[0].structure_energies, np.ndarray)
        assert isinstance(self.outcomes[0].dofs_energies, np.ndarray)
        assert isinstance(self.outcomes[0].dim_indices, np.ndarray)
        assert isinstance(self.outcomes[0].structures, list)

    def _populate_variables(self, kindanalysis):
        """
        This method creates the attributes used fot the analysis extracting the
        values from the sith_objects.
        """
        self.all_dft_energies = []
        self.all_edm_energies = [] # energy distribution method
        for ed_obj in self.outcomes:
            self.all_dft_energies.append(ed_obj.structures_scf_energies)
            self.all_edm_energies.append(ed_obj.structure_energies)

        self.all_dft_energies = np.array(self.all_dft_energies)
        self.all_edm_energies = np.array(self.all_edm_energies)

        return self

    def plot_le(self, a_names, aminos=3, ax: plt.Axes = None, sp=None,
                lw=1, ms=1, **kwargs):
        """
        plots the 
        """
        if 'ax_pref' in kwargs:
            setter = kwargs['ax_pref']
            del kwargs['ax_pref']
        else:
            setter = {}

        if sp is None:
            sp = StandardPlotter(**kwargs)
        if ax is None:
            ax = sp.ax[0]
        sp.axis_setter(ax=ax,
                        xlabel=f'Distance({", ".join(a_names)})[\u212B]',
                        ylabel='Energy[Ha]',
                        **setter)
        ls = []
        es = []
        for an in self.analysis:
            l, e = an.le_dof_amino(a_names, aminos)
            ls.append(l)
            es.append(e)
            sp.plot_data(l, e, ax=ax, lw=lw, markersize=ms)
        return ax, ls, es
    
    def le_all(self, a_names, aminos):
        self.ls = []
        self.es = []
        for an in self.analysis:
            l, e = an.le_dof_amino(a_names, aminos)
            self.ls.append(l)
            self.es.append(e)
        


    def plot_DFT_ener(self, ax: plt.Axes = None, sp=None,
                      lw=1, ms=1, **kwargs):
        """
        plots the 
        """
        setter = kwargs['ax_pref']
        del kwargs['ax_pref']

        if sp is None:
            sp = StandardPlotter(**kwargs)
        if ax is None:
            ax = sp.ax[0]
        sp.axis_setter(ax=ax,
                       xlabel=f'Stretched Structure',
                       ylabel='DFT Energy[Ha]',
                       **setter)
        
        xs = []
        ys = []
        for i, sith in enumerate(self.outcomes):
            y = sith.structures_scf_energies
            index1 = self.pep_infos[i].amino_info[1]['CH3'] - 1
            index2 = self.pep_infos[i].amino_info[5]['CH3'] - 1
            x = []
            for struc in sith.structures:
                dist = struc.atoms.get_distance(index1, index2)
                x.append(dist)
            sp.plot_data(x, y, ax=ax, lw=lw, markersize=ms)
            xs.append(x)
            ys.append(y)

        return ax, xs, ys
    
    def amino_freq(self):    
        self.names = []
        for sith in self.outcomes:
            self.names.append(list(sith.name))
        self.names = np.array(self.names)

        return self.names
