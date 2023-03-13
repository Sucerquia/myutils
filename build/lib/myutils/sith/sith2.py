import numpy as np
from ase import units
import glob


class void:
    pass


class Sith2:
    def __init__(self, sith, index_file=None, forces_files=None):
        self.sith = sith
        self.read_forces_and_indexes(index_file, forces_files)
        self.all_forces = self.extract_forces()
        self.delta_rics = self.extract_changes()
        self._reference = self._deformed[0]

    def read_forces_and_indexes(self, index_file=None, forces_files=None):
        if forces_files is None:
            forces_files = glob.glob('./forces/*force*.dat')

        if index_file is None:
            index_file = './forces/indexes.dat'

        self.atomic_symbols = np.loadtxt(index_file, usecols=1, dtype=str)
        assert (self.atomic_symbols == self.sith._reference.atoms.get_chemical_symbols()).all()

        forces_files.sort()
        self._deformed = [void() for i in forces_files]

        print(" --- reading forces ----")
        for i, file in enumerate(forces_files):
            print(i, file)
            self._deformed[i].internal_forces = np.loadtxt(file, usecols=2)
            indexes = np.loadtxt(file, usecols=1, dtype=str)
            self._deformed[i].dimIndices = [eval(ric_indexes) for ric_indexes in indexes]
            self._deformed[i].ric = np.loadtxt(file, usecols=3)
        self.dims=np.zeros((4), dtype=int)
        for index in self._deformed[0].dimIndices:
            self.dims[0] += 1
            self.dims[len(index) - 1] += 1
        for defo in self._deformed:
            defo.dims = self.dims

        self.indexes = self._deformed[0].dimIndices

        return self.indexes, self._deformed

    def extract_changes(self):
        rics = list()
        for defo in self._deformed:
            ric = defo.ric
            ric[:self.dims[1]] = ric[:self.dims[1]] / units.Bohr
            ric[self.dims[1]:] = ric[self.dims[1]:] * np.pi/180 
            rics.append(defo.ric)

        rics = np.array(rics)
        delta_rics = rics - np.insert(rics[:-1], 0, rics[0], axis=0)

        with np.nditer(delta_rics, op_flags=['readwrite']) as dqit:
            for dq in dqit:
                dq[...] = dq - 2*np.pi if dq > np.pi \
                    else (dq + 2*np.pi if dq < -np.pi else dq)

        return delta_rics

    def extract_forces(self):
        all_forces = list()
        for defo in self._deformed:
            all_forces.append(defo.internal_forces)

        return np.array(all_forces)

    def analysis(self):
        all_values = - self.all_forces * self.delta_rics
        self.energies = np.cumsum(all_values, axis=0)
        self.total_ener = np.sum(self.energies, axis=1)

        return self.energies, self.total_ener

