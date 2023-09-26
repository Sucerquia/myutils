from SITH.SITH import SITH
import numpy as np
from ase import units


class Jedi(SITH):
    """
    This tool takes the previous created sith and keeps some things in order
    like:
    - same units always (Ha, A, rad)
    - all objects (all_rics, all_hessians...)
    """
    def __init__(self, drelaxed, dstreched, atoms2kill=None,
                 rem_first_def=0, rem_last_def=0, from_last_minima=True):
        SITH.__init__(self, drelaxed, dstreched)
        if atoms2kill is not None:
            self.setKillAtoms(atoms2kill)
        self.extractData()
        self.dims = self._deformed[0].dims
        
        self.all_rics = self.rics() 
        self.all_rics[:, :self.dims[1]] *= units.Bohr
        self.all_hessians = np.array([defo.hessian for defo in self._deformed])
        self.all_hessians[:, :, :self.dims[1]] /= units.Bohr
        self.all_hessians[:, :self.dims[1]] /= units.Bohr
        self.n_deformed = len(self._deformed)
        self.scf_energy = np.array([defo.energy - self._deformed[0].energy
                                    for defo in self._deformed])
        if from_last_minima:
            try:
                dif_ener = self.scf_energy[1:] - self.scf_energy[:-1]
                rem_first_def = (np.where(dif_ener < 0)[0] + 1)[-1]
            except IndexError:
                pass

        # the next step creates
        # - self.deltaQ
        # - self.energies
        # - self.deformationEnergy
        self.rem_first_last(rem_first_def, rem_last_def)

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
            ric[defo.dims[1]:] = ric[defo.dims[1]:]
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
        delta_rics = self.all_rics - self.all_rics[0]

        condition = delta_rics[:, self.dims[1]:] > np.pi
        delta_rics[:, self.dims[1]:][condition] -= 2 * np.pi
        condition = delta_rics[:, self.dims[1]:] < -np.pi
        delta_rics[:, self.dims[1]:][condition] += 2 * np.pi

        return delta_rics

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
        ini_index = rem_first_def
        last_index = self.n_deformed - rem_last_def
        self._deformed = self._deformed[ini_index: last_index]
        self.scf_energy = self.scf_energy[ini_index: last_index] - \
            self.scf_energy[ini_index]
        self.all_rics = self.all_rics[ini_index: last_index]
        self.n_deformed -= rem_first_def + rem_last_def
        self.all_hessians = self.all_hessians[ini_index: last_index]
        self.deltaQ = self.extract_changes()
        self.energies = self.energyAnalysis()
        self.deformationEnergy = np.sum(self.energies, axis=1)

        return self._deformed

    def energyAnalysis(self):
        """
        Performs the SITH energy analysis, populates energies,
        deformationEnergy, and pEnergies.

        Notes
        -----
        Consists of the dot multiplication of the deformation vectors and the
        Hessian matrix (analytical gradient of the harmonic potential energy
        surface) to produce both the total calculated change in energy between
        the reference structure and each deformed structure
        (SITH.deformationEnergy) as well as the subdivision of that energy into
        each DOF (SITH.energies).
        """
        print("Performing energy analysis...")
        if self.deltaQ is None or self.q0 is None or self.qF is None:
            raise Exception(
                "Populate Q has not been executed so necessary data for "
                "analysis is lacking. This is likely due to not calling "
                "extractData().")
        self.energies = 0.5 * np.matmul(self.deltaQ, self.all_hessians[0]) * \
            self.deltaQ

        print("Execute Order 67. Successful energy analysis completed.")
        return self.energies

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
