import numpy as np
from ase.io import read
from myutils.sith.visualization import MoleculeViewer


class info:
    def __init__(self, pdb_file):
        self.atoms = read(pdb_file)
        indexes_aminos = self.atoms.arrays['residuenumbers']
        names_aminos = self.atoms.arrays['residuenames']
        indexes_atoms = np.arange(len(self.atoms)) + 1
        self.name_atoms_raw = self.atoms.arrays['atomtypes']

        self.atom_names = {}
        self.amino_name = {}
        self.atom_indexes = {}
        self.amino_info = {}

        for i in range(len(indexes_aminos)):
            self.amino_name[indexes_aminos[i]] = names_aminos[i]
            self.atom_names[indexes_aminos[i]] = []
            self.atom_indexes[indexes_aminos[i]] = []
            self.amino_info[indexes_aminos[i]] = {}
        for i in range(len(indexes_aminos)):
            index = indexes_aminos[i]
            self.atom_names[index].append(self.name_atoms_raw[i])
            self.atom_indexes[index].append(indexes_atoms[i])
            self.amino_info[index][self.name_atoms_raw[i]] = indexes_atoms[i]

    def endo_exo_proline(self, sith):
        """
        Endo, exo states must correspond with the first angle equal to zero.
        """
        prolines = np.where(np.array(
                            list(self.amino_name.values())) == 'PRO')[0]
        angles = np.array([[0, 0]])
        for i in prolines:
            isolated_proline = self.amino_info[i + 1]
            for conf in sith._deformed:
                try:
                    index = conf.dimIndices.index((isolated_proline['CB'],
                                                   isolated_proline['CA'],
                                                   isolated_proline['N'],
                                                   isolated_proline['CD']))
                    angle1 = conf.ric[index]*180/np.pi
                except ValueError:
                    index = conf.dimIndices.index((isolated_proline['CD'],
                                                   isolated_proline['N'],
                                                   isolated_proline['CA'],
                                                   isolated_proline['CB']))
                    angle1 = -conf.ric[index]*180/np.pi
                try:
                    index = conf.dimIndices.index((isolated_proline['N'],
                                                   isolated_proline['CA'],
                                                   isolated_proline['CB'],
                                                   isolated_proline['CG']))
                    angle2 = conf.ric[index]*180/np.pi
                except ValueError:
                    index = conf.dimIndices.index((isolated_proline['CG'],
                                                   isolated_proline['CB'],
                                                   isolated_proline['CA'],
                                                   isolated_proline['N']))
                    angle2 = -conf.ric[index]*180/np.pi
                angles = np.append(angles, [[angle1, angle2]], axis=0)

        return np.delete(angles, 0, 0)

    def xis_angles(self, sith):
        """
        Endo, exo states must correspond with the first angle equal to zero.
        """
        prolines = np.where(np.array(
                            list(self.amino_name.values())) == 'PRO')[0]
        angles = np.array([[0, 0]])
        for i in prolines:
            isolated_proline = self.amino_info[i + 1]
            for conf in sith._deformed:
                try:
                    index = conf.dimIndices.index((isolated_proline['N'],
                                                   isolated_proline['CD'],
                                                   isolated_proline['CG'],
                                                   isolated_proline['CB']))
                    angle1 = -conf.ric[index]*180/np.pi
                except ValueError:
                    index = conf.dimIndices.index((isolated_proline['CB'],
                                                   isolated_proline['CG'],
                                                   isolated_proline['CD'],
                                                   isolated_proline['N']))
                    angle1 = conf.ric[index]*180/np.pi
                try:
                    index = conf.dimIndices.index((isolated_proline['N'],
                                                   isolated_proline['CA'],
                                                   isolated_proline['CB'],
                                                   isolated_proline['CG']))
                    angle2 = conf.ric[index]*180/np.pi
                except ValueError:
                    index = conf.dimIndices.index((isolated_proline['CG'],
                                                   isolated_proline['CB'],
                                                   isolated_proline['CA'],
                                                   isolated_proline['N']))
                    angle2 = -conf.ric[index]*180/np.pi
                angles = np.append(angles,
                                   [[(angle1+angle2)/(2*np.cos(4*np.pi/5)),
                                     (angle1-angle2)/(2*np.cos(4*np.pi/5))]],
                                   axis=0)

        return np.delete(angles, 0, 0)

    def compute_dihedrals(self, atom1index, atom2index, atom3index,
                          atom4index):
        # axis of rotation of the dihedral angle
        axis = (self.atoms[atom3index-1].position -
                self.atoms[atom2index-1].position)
        # vector bewteen atoms 1 and 2
        axis1 = (self.atoms[atom1index-1].position -
                 self.atoms[atom2index-1].position)
        # vector bewteen atoms 3 and 4
        axis2 = (self.atoms[atom4index-1].position -
                 self.atoms[atom3index-1].position)

        # initial and final sides of the angle in the same plane
        a = axis1 - axis * (np.dot(axis, axis1)/np.dot(axis, axis))
        b = axis2 - axis * (np.dot(axis, axis2)/np.dot(axis, axis))

        len_axis = np.linalg.norm(axis)
        lena = np.linalg.norm(a)
        lenb = np.linalg.norm(b)
        # value of the angle
        theta = np.arccos(np.dot(a, b)/(lena * lenb))

        # sign
        cross = axis + np.cross(a, b)
        len_cross = np.linalg.norm(cross)
        sign = (len_cross - len_axis)/abs(len_cross - len_axis)

        return sign * theta

    def build_proline_state(self, aminoacid, state):
        """
        Parameters
        ==========
        aminoacid: int
            number of the aminoacid to modify. It has to correspond
            with a proline aminoacid.

        state: str
            'endo' or 'exo'. States of the proline
        """
        # test for proline
        name_atoms = list(self.amino_info[aminoacid].keys())
        assert len(name_atoms) == 14, "wrong amino acid, check that you are" +\
            " trying to define a proline and not another one."

        proline_atoms = ['N', 'CD', 'HD1', 'HD2', 'CG', 'HG1', 'HG2', 'CB',
                         'HB1', 'HB2', 'CA', 'HA', 'C', 'O']
        for atom in proline_atoms:
            assert atom in name_atoms

        # defining ring atoms indexes:
        ca = self.amino_info[aminoacid]['CA'] - 1
        cb = self.amino_info[aminoacid]['CB'] - 1
        cg = self.amino_info[aminoacid]['CG'] - 1
        hg1 = self.amino_info[aminoacid]['HG1'] - 1
        hg2 = self.amino_info[aminoacid]['HG2'] - 1
        cd = self.amino_info[aminoacid]['CD'] - 1
        hd1 = self.amino_info[aminoacid]['HD1'] - 1
        hd2 = self.amino_info[aminoacid]['HD2'] - 1
        n = self.amino_info[aminoacid]['N'] - 1

        # align Ca N  CB in the same plane with Ca N in the x axis
        v = MoleculeViewer(self.atoms, alignment=[ca, n, cb])
        self.atoms = v.atoms
        # move Cd to the xy plane
        vec = self.atoms[cd].position
        trans = v.align_plane(vec)
        v.apply_trans(self.atoms, trans, indexes=[cd, hd1, hd2])

        # move Cg to the xy plane
        vec = self.atoms[cg].position
        trans = v.align_plane(vec)

        v.apply_trans(self.atoms, trans, indexes=[cg, hg1, hg2])

        # align cb, cd in the x axis, cg in the xy plane
        v.xy_alignment(self.atoms, cb, cd, cg)

        # endo-exo state
        if state == 'endo':
            bending = -30 * np.pi/180
        elif state == 'exo':
            bending = 30 * np.pi/180
        else:
            raise ValueError("The state of proline has to be either endo" +
                             " or exo")

        v.apply_trans(self.atoms, v.rot_x(bending), indexes=[cg, hg1, hg2])

        return self.atoms
