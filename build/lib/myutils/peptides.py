import numpy as np


class info:
    def __init__(self, pdb_file, withx=0):
        pdb_header = ['END',
                      'MODEL',
                      'CRYST1',
                      'REMARK',
                      'TITLE',
                      'TER']
        indexes_aminos = np.loadtxt(pdb_file, comments=pdb_header,
                                    usecols=5+withx, dtype=int)
        names_aminos = np.loadtxt(pdb_file, comments=pdb_header,
                                  usecols=3, dtype=str)
        indexes_atoms = np.loadtxt(pdb_file, comments=pdb_header,
                                   usecols=1, dtype=int)
        self.name_atoms_raw = np.loadtxt(pdb_file, comments=pdb_header,
                                         usecols=2, dtype=str)
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
