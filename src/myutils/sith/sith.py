import numpy as np
from ase import units
import glob
from myutils.miscellaneous import output_terminal
from ase.io import read
from scipy.integrate import simpson


class Geometry:
    def __init__(self, xyz_file, force_file):
        self.name = force_file
        self.atoms = read(xyz_file)
        self.nAtoms = len(self.atoms)
        indexes = np.loadtxt(force_file, usecols=1,
                             dtype=str)
        self.dimIndices = [eval(ric_indexes) for ric_indexes in indexes]
        self.ric = np.loadtxt(force_file, usecols=3)

        self.dims = np.zeros((4), dtype=int)
        for index in self.dimIndices:
            self.dims[0] += 1
            self.dims[len(index) - 1] += 1

        self.internal_forces = np.loadtxt(force_file, usecols=2)
        self.internal_forces[:self.dims[1]] = \
            self.internal_forces[:self.dims[1]] / units.Bohr

        with open(force_file) as ff:
            header = ff.readline()
        self.energy = float(header.replace('\n', '').split(' ')[-1])

    def _killDOFs(self, dofis):
        """
        Takes a list of indices of degrees of freedom to remove, Removes DOFs
        from ric, dimIndices, and internal_forces, updates dims
        """
        self.ric = np.delete(self.ric, dofis)
        # counter of the number of DOF removed arranged in
        # types (lenght, angle, dihedral)
        tdofremoved = [0, 0, 0]

        for index in sorted(dofis, reverse=True):
            tdofremoved[len(self.dimIndices[index]) - 2] += 1
            del self.dimIndices[index]

        self.dims[0] -= len(dofis)
        self.dims[1] -= tdofremoved[0]
        self.dims[2] -= tdofremoved[1]
        self.dims[3] -= tdofremoved[2]
        # change for forces
        if(self.internal_forces is not None):
            self.internal_forces = np.delete(self.internal_forces, dofis)
        return self.ric, self.dims, self.internal_forces


class Sith:
    def __init__(self, forces_xyz_files=None, master_directory='./',
                 killAtoms=None, killDOFs=None, killElements=None,
                 rem_first_def=0, rem_last_def=0, method=2):

        # General variabales
        # Define files
        self.setting_force_xyz_files(forces_xyz_files, master_directory)
    def setting_force_xyz_files(self, forces_xyz_files=None,
                                master_directory='./'):
        """
        Creates the list of forces and xyz files.

        Parameters
        ==========
        master_directory: str
            path to the directory containing all the necessary files.
        forces_xyz_files: list[2 str] (optional)
            path to the forces files and the xyz files.
            Default=['master_directory/*force*.dat,
                     'master_directory/*force*.xyz]

        Return
        ======
        (tuple) [2list, #deformed] List of forces files (first element) and
        list of xyz files (second element).
        """
        if forces_xyz_files is None:
            forces_xyz_files = [None, None]
        if killAtoms is None:
            killAtoms = []
        if killDOFs is None:
            killDOFs = []
        if killElements is None:
            killElements = []

        # files path
        self.forces_files = forces_xyz_files[0]
        self.xyz_files = forces_xyz_files[1]
        self.master_directory = master_directory

        # default set-up
        if self.forces_files is None:
            self.forces_files = glob.glob(self.master_directory +
                                          '/*force*.dat')
        if self.xyz_files is None:
            self.xyz_files = glob.glob(master_directory + '/*force*.xyz')
        assert len(self.forces_files) == len(self.xyz_files), "Different " + \
            "number of forces and xyz files."
        # # Create forces files
        if (len(self.forces_files) == 0):
            if (len(glob.glob(self.master_directory+'/*force*.log')) == 0):
                raise OSError(f"{self.master_directory} does not exist or " +
                              "does not contain *force*.log files the " +
                              "forces_xyz_files is empty.")
            else:
                self.create_files()

        # # sorth forces and xyz
        self.forces_files.sort()
        self.xyz_files.sort()

        # create Geometries shape=(n_def, 1)
        self._deformed = [Geometry(i, j)
                          for i, j in zip(self.xyz_files, self.forces_files)]
        # number of deformed configs
        self.n_deformed = len(self._deformed)

        # debug: test the DOFs in all force files
        self.check_dofs()

        # Create matrices
        # # DFT energies for each configuration shape=(n_def, 1)
        self.deformationEnergy = np.array([defo.energy
                                           for defo in self._deformed])
        # # ric values matrix shape=(n_def, n_dofs)
        self.qF = np.array([defo.ric for defo in self._deformed])
        # # matrix changes shape=(n_def, n_dofs)
        self.all_rics = self.rics()
        self.deltaQ = self.extract_changes()
        # # all forces shape=(n_def, n_dofs)
        self.all_forces = np.array([defo.internal_forces
                                    for defo in self._deformed])

        # # energies per DOF shape=(n_def, n_dofs)
        # # and computed energy shape=(n_def, 1)
        if self.method == 1:
            self.energies, self.configs_ener = self.analysis()
        elif self.method == 2:
            self.energies, self.configs_ener = self.analysis_classical()
        else:
            raise ValueError("Non-recognized method")

        self.killer(killAtoms, killDOFs, killElements)
        self.rem_first_last(rem_first_def, rem_last_def)

    def create_files(self):
        """"
        Uses \"myutils extract_forces\" to get the forces, coordinates and
        coordinates in xyz format.
        """
        get_forces_exec = output_terminal("myutils extract_forces")
        output_terminal(get_forces_exec.replace("\n", "") +
                        f" -d {self.master_directory}")
        self.forces_files = glob.glob(self.master_directory + '/*force*.dat')
        self.xyz_files = glob.glob(self.master_directory + '/*force*.xyz')

        return self.forces_files, self.xyz_files

    def check_dofs(self):
        referece = self._deformed[0].dimIndices
        for defo in self._deformed[1:]:
            to_compare = defo.dimIndices
            assert len(to_compare) == len(referece), \
                "The number of DOFs in the first force file and in " + \
                f"{defo.name} is different."
            for i in range(len(referece)):
                assert to_compare[i] == referece[i] or \
                    to_compare[i] == referece[i][::-1], \
                    f"The DOF in the first force file is {referece[i]} and" +\
                    f" the DOF in {defo.name} is {to_compare[i]}"
        return True

    def rics(self):
        rics = list()
        for defo in self._deformed:
            ric = defo.ric
            ric[:defo.dims[1]] = ric[:defo.dims[1]]
            ric[defo.dims[1]:] = ric[defo.dims[1]:] * np.pi/180
            rics.append(defo.ric)
        return np.array(rics)

    def extract_changes(self):
        delta_rics = self.all_rics - np.insert(self.all_rics[:-1], 0,
                                               self.all_rics[0],
                                               axis=0)

        with np.nditer(delta_rics, op_flags=['readwrite']) as dqit:
            for dq in dqit:
                dq[...] = dq - 2*np.pi if dq > np.pi \
                    else (dq + 2*np.pi if dq < -np.pi else dq)

        return delta_rics

    def analysis(self):
        all_ener = [np.zeros(len(self.all_forces[0]))]
        # next loop is a nasty cummulative integration. Maybe it could
        # be improved
        for i in range(1, len(self.all_forces) + 1):
            ener_def = -simpson(self.all_forces[:i + 1],
                                x=self.all_rics[:i + 1],
                                axis=0)
            all_ener.append(ener_def)
        all_ener = np.array(all_ener)
        total_ener = np.sum(all_ener, axis=1)

        return all_ener, total_ener

    def analysis_classical(self):
        all_values = - self.all_forces * self.deltaQ
        energies = np.cumsum(all_values, axis=0)
        total_ener = np.sum(energies, axis=1)

        return energies, total_ener

    def compareEnergies(self):
        """
        Takes in SITH object sith, Returns Tuple of expected stress energy,
        stress energy error, and %Error

        Notes
        -----
        Expected Stress Energy: Total E deformed structure from input .fchk -
        total E reference structure from input .fchk
        Stress Energy Error: calculated stress energy - Expected Stress Energy
        %Error: Stress Energy Error / Expected Stress Energy
        """
        obtainedDE = self.configs_ener
        expectedDE = self.deformationEnergy - self.deformationEnergy[0]

        errorDE = obtainedDE - expectedDE
        pErrorDE = (errorDE / expectedDE) * 100
        pErrorDE[0] = 0
        pErrorDE[np.logical_not(np.isfinite(pErrorDE))] = 200

        return np.array([obtainedDE, expectedDE, errorDE, pErrorDE])

# --------------------------- killer ------------------------------------------
    def killer(self, killAtoms=[], killDOFs=[], killElements=[]):
        """
        killAtoms:
            list of indexes of atoms to be killed
        killDOFs:
            list of tuples with the DOFs to be killed
        killElements:
            list of strings with the elements to be killed. So, if you want to
            remove all hydrogens and carbons, use killElements=['H', 'C'].
        """

        self.dims_to_kill = killDOFs
        self.atoms_to_kill = killAtoms

        # concatenate elements in atoms to be killed
        molecule = np.array(self._deformed[0].atoms.get_chemical_symbols())

        for element in killElements:
            indexes_element = np.where(molecule == element)[0] + 1
            self.atoms_to_kill.extend(indexes_element)

        # concatenate atoms in DOFs to be killed
        for atom in self.atoms_to_kill:
            self.dims_to_kill.extend(
                [dim for dim in self._deformed[0].dimIndices if atom in dim])
        # remove repetitions
        self.dims_to_kill = list(dict.fromkeys(self.dims_to_kill))

        # remove DOFs
        self.__killDOFs(self.dims_to_kill)

        if self.method == 1:
            self.energies, self.configs_ener = self.analysis()
        elif self.method == 2:
            self.energies, self.configs_ener = self.analysis_classical()
        else:
            raise ValueError("Non-recognized method")

        return self.dims_to_kill

    def __killDOFs(self, dofs):
        rIndices = list()
        for dof in dofs:
            rIndices.extend([i for i in range(self._deformed[0].dims[0])
                            if self._deformed[0].dimIndices[i] == dof])

        # kill DOFs in Geometries
        for defo in self._deformed:
            defo._killDOFs(rIndices)

        # kill DOFs in sith
        self.qF = np.delete(self.qF,  rIndices, axis=1)
        self.all_rics = np.delete(self.all_rics,  rIndices, axis=1)
        self.deltaQ = np.delete(self.deltaQ,  rIndices, axis=1)
        self.all_forces = np.delete(self.all_forces,  rIndices, axis=1)

        return rIndices

# ---------------------- remove deformations ----------------------------------
    def rem_first_last(self, rem_first_def=0, rem_last_def=0):
        self._deformed = self._deformed[rem_first_def:self.n_deformed -
                                        rem_last_def]
        self.qF = self.qF[rem_first_def:self.n_deformed -
                          rem_last_def]
        self.deltaQ = self.deltaQ[rem_first_def: self.n_deformed -
                                  rem_last_def]
        self.all_forces = self.all_forces[rem_first_def: self.n_deformed -
                                          rem_last_def]
        self.energies = self.energies[rem_first_def: self.n_deformed -
                                      rem_last_def]
        self.configs_ener = self.configs_ener[rem_first_def: self.n_deformed -
                                              rem_last_def]
        self.deformationEnergy = self.deformationEnergy[rem_first_def:
                                                        self.n_deformed -
                                                        rem_last_def]
        self.all_rics = self.all_rics[rem_first_def:
                                      self.n_deformed -
                                      rem_last_def]

        return self._deformed
