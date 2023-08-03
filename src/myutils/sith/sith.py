import numpy as np
from ase import units
import glob
from myutils.miscellaneous import output_terminal
from ase.io import read
from scipy.integrate import simpson


class Geometry:
    def __init__(self, xyz_file, force_file):
        """
        A class that contains the information of one stretched configuration
        for the sith analysis.

        Parameters
        ==========
        xyz_file: str
            config file of the stretched config.
        force_file: str
            forces data file of the stretched config.
        """
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
        Takes a list of indices of degrees of freedom to be removed and removes
        them from ric, dimIndices, and internal_forces, and updates dims.

        Parameters
        ==========
        dofis: list[int]
            list of indexes to be removed.

        Return
        ======
        (tuple) [DIFs, dims, intertanl_forces] removing the desired DOFs.
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
                 rem_first_def=0, rem_last_def=0, integration_method=1,
                 name=None):
        """
        A class to calculate & house SITH analysis data.

        Parameters
        ==========
        master_directory: str
            path to the directory containing all the forces.dat files
        forces_xyz_files: list[2 str]
            path to the forces files and the xyz files.
            Default =['master_directory/*force*.dat,
                     'master_directory/*force*.xyz]
        killAtoms: list[ints]
            List of indices of atoms to be removed of the analysis. Typically
            light atoms. Default=None
        killDOFs: list[ int ]
            indexes of degrees of freedom to be removed of the analysis.
            Default=None
        killElements: list[str]
            chemical symbols of atoms to be removed of the analysis.
            Default=None
        rem_first_def: int
            number of first deformations to be removed. Used when the system
            suffer large changes at the begining. Default=0
        rem_last_def: int
            number of last deformations to be removed. Used when a rupture is
            produced. Default=0
        integration_method: int
            Index of numertical integration method according to the list
            [self.rectangle_integration, self.trapezoid_integration,
            self.simpson_integration]. Default=1 (namely, trapezoid)
        name: str
            name of the molecule to be analyze. This is completely arbitrary
            and up to the user. Default=master_directory

        Note: this code is done such that the deformation is sorted in
        alphabetic order.
        """
        if name is None:
            name = master_directory
        self.name = name

        # Define files
        self.setting_force_xyz_files(forces_xyz_files, master_directory)

        # create Geometries shape=(n_def, 1)
        self._deformed = [Geometry(i, j)
                          for i, j in zip(self.xyz_files, self.forces_files)]
        # number of deformed configs
        self.n_deformed = len(self._deformed)
        self.dims = self._deformed[0].dims

        # debug: test the DOFs in all force files
        self._check_dofs()

        # Create matrices
        # # DFT energies for each configuration shape=(n_def, 1)
        self.scf_energy = np.array([defo.energy - self._deformed[0].energy
                                    for defo in self._deformed])
        # # ric values matrix shape=(n_def, n_dofs)
        self.qF = np.array([defo.ric for defo in self._deformed])
        self.all_rics = self.rics()
        # # matrix changes shape=(n_def, n_dofs)
        self.deltaQ = self.extract_changes()
        # # all forces shape=(n_def, n_dofs)
        self.all_forces = np.array([defo.internal_forces
                                    for defo in self._deformed])

        # Numerical integration
        # # energies per DOF shape=(n_def, n_dofs)
        # # and computed energy shape=(n_def, 1)
        implemented_methods = [self.rectangle_integration,
                               self.trapezoid_integration,
                               self.simpson_integration]
        self.integration_method = implemented_methods[integration_method]
        self.energies, self.deformationEnergy = self.integration_method()

        # remove atoms, dofs and last or first configurations
        if killAtoms is None:
            killAtoms = []
        if killDOFs is None:
            killDOFs = []
        if killElements is None:
            killElements = []

        self.killer(killAtoms, killDOFs, killElements)
        self.rem_first_last(rem_first_def, rem_last_def)

    def setting_force_xyz_files(self, forces_xyz_files=None,
                                master_directory='./'):
        """
        Creates the list of forces and xyz files.

        Parameters
        ==========
        master_directory: str
            path to the directory containing all the necessary files.
            Default='./'
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

        # Create forces files
        if (len(self.forces_files) == 0):
            if (len(glob.glob(self.master_directory + '/*force*.log')) == 0):
                raise OSError(f"{self.master_directory} does not exist or " +
                              "does not contain *force*.log")
            else:
                self.create_files()

        # # sorth forces and xyz
        self.forces_files.sort()
        self.xyz_files.sort()

        return self.forces_files, self.xyz_files

    def create_files(self):
        """"
        Uses \"myutils extract_forces\" to get the forces, coordinates and
        coordinates in xyz format from the master directory.

        Return
        ======
        (tuple) [2list, #deformed] List of forces files (first element) and
        list of xyz files (second element).
        """
        get_forces_exec = output_terminal("myutils extract_forces")
        output_terminal(get_forces_exec.replace("\n", "") +
                        f" -d {self.master_directory}")
        self.forces_files = glob.glob(self.master_directory + '/*force*.dat')
        self.xyz_files = glob.glob(self.master_directory + '/*force*.xyz')

        self.forces_files.sort()
        self.xyz_files.sort()

        return self.forces_files, self.xyz_files

    def _check_dofs(self):
        """
        Checks the DOFs if all forces files are the same.

        Return
        ======
        (bool) [True] In case it is successful. If it fails, it raises an
        assert error.
        """
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
        """
        Extract and concatenate the DOFs values in a matrix. Angles are given
        in radians.

        Return
        ======
        (np.array) [#def, #dofs] matrix of DOFs values per deformed
        configuration.
        """
        rics = list()
        for defo in self._deformed:
            ric = defo.ric.copy()
            ric[defo.dims[1]:] = ric[defo.dims[1]:] * np.pi / 180
            rics.append(ric)
        return np.array(rics)

    def extract_changes(self):
        """
        Extract and concatenate the DOFs changes in a matrix. Angles are given
        in radians.

        Return
        ======
        (np.array) [#def, #dofs] matrix of DOFs changes per deformed
        configuration.
        """
        delta_rics = self.all_rics - np.insert(self.all_rics[:-1], 0,
                                               self.all_rics[0],
                                               axis=0)

        delta_rics[:, self.dims[1]:][delta_rics[:, self.dims[1]:]
                                     > np.pi] -= 2 * np.pi
        delta_rics[:, self.dims[1]:][delta_rics[:, self.dims[1]:]
                                     < -np.pi] += 2 * np.pi

        return delta_rics

    def rectangle_integration(self):
        """
        Numerical integration using rectangle rule algorithm. Method 0 in this
        class (see Sith parameters).

        Return
        ======
        (tuple) [energies, total_ener] energies computed by SITH
        method.
        """
        all_values = - self.all_forces * self.deltaQ
        energies = np.cumsum(all_values, axis=0)
        total_ener = np.sum(energies, axis=1)

        return energies, total_ener

    def trapezoid_integration(self):
        """
        Numerical integration using trapezoid rule algorithm. Method 1 in this
        class (see Sith parameters).

        Return
        ======
        (tuple) [energies, total_ener] energies computed by SITH
        method.
        """
        # energy for the optimized config must be the reference

        added_forces = (self.all_forces[1:] + self.all_forces[:-1]) / 2
        all_values = added_forces * self.deltaQ[1:]
        all_values = np.insert(all_values, 0, np.zeros(self.dims[0]), axis=0)
        energies = -np.cumsum(all_values, axis=0)
        total_ener = np.sum(energies, axis=1)

        return energies, total_ener

    def simpson_integration(self):
        """
        Numerical integration using simpson algorithm. Method 2 in this class
        (see Sith parameters).

        Return
        ======
        (tuple) [energies, total_ener] energies computed by SITH
        method.
        """
        rics = self.all_rics.copy().T

        for angle in rics[self.dims[1]:]:
            for i in range(self.n_deformed - 1):
                while angle[i + 1] - angle[i] > np.pi:
                    angle[i + 1:] -= 2 * np.pi
                while angle[i + 1] - angle[i] < -np.pi:
                    angle[i + 1:] += 2 * np.pi
        rics = rics.T

        # first array counts the  energy in the dofs for the optimized
        # configuration. that's why it is zero
        all_ener = np.array([[0] * self.dims[0]])
        # next loop is a 'nasty' cummulative integration. Maybe it could
        # be improved
        for i in range(1, self.n_deformed):
            ener_def = -simpson(self.all_forces[: i + 1],
                                x=rics[: i + 1],
                                axis=0)
            all_ener = np.append(all_ener, np.array([ener_def]), axis=0)
        total_ener = np.sum(all_ener, axis=1)
        return all_ener, total_ener

    def compareEnergies(self):
        """
        computes the difference between the expected and the obtained change on
        energy.

        Return
        ======
        (np.array)[energy, expected energy, error, percentage of error]
        """
        obtainedDE = self.deformationEnergy
        expectedDE = self.scf_energy

        errorDE = obtainedDE - expectedDE
        expectedDE[abs(expectedDE) < 1e-15] = 1e-15
        pErrorDE = (errorDE / expectedDE) * 100
        pErrorDE[0] = 0
        pErrorDE[np.logical_not(np.isfinite(pErrorDE))] = 200

        return np.array([obtainedDE, expectedDE, errorDE, pErrorDE])

# --------------------------- killer ------------------------------------------
    def killer(self, killAtoms=[], killDOFs=[], killElements=[]):
        """
        Removes all DOFs who contain atoms, DOFs and elements to be removed.

        Parameters
        ==========
        killAtoms:
            list of indexes of atoms to be killed
        killDOFs:
            list of tuples with the DOFs to be killed
        killElements:
            list of strings with the elements to be killed. So, if you want to
            remove all hydrogens and carbons, use killElements=['H', 'C'].

        Return
        ======
        (list) [int] indexes of the removed DOFs.
        """

        self.dims_to_kill = [self._deformed[0].dimIndices[i] for i in killDOFs]
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
        self.qF = np.delete(self.qF, rIndices, axis=1)
        self.all_rics = np.delete(self.all_rics, rIndices, axis=1)
        self.deltaQ = np.delete(self.deltaQ, rIndices, axis=1)
        self.all_forces = np.delete(self.all_forces, rIndices, axis=1)

        return rIndices

# ---------------------- remove deformations ----------------------------------
    def rem_first_last(self, rem_first_def=0, rem_last_def=0):
        """
        Removes first and last deformation configs and data from all the
        attributes of the Sith object.

        Parameters
        ==========
        rem_first_def: int
            number of configuration to remove in the first stretched
            configuration.
        rem_last_def: int
            number of configuration to remove in the last stretching
            configuration

        Return
        ======
        (list) Deformed Geometry objects. self._deformed.
        """
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
        self.deformationEnergy = self.deformationEnergy[rem_first_def:
                                                        self.n_deformed -
                                                        rem_last_def]
        self.scf_energy = self.scf_energy[rem_first_def:
                                          self.n_deformed -
                                          rem_last_def]
        self.all_rics = self.all_rics[rem_first_def:
                                      self.n_deformed -
                                      rem_last_def]
        self.n_deformed -= rem_first_def + rem_last_def

        return self._deformed
