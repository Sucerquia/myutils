import numpy as np
from ase.calculators.gaussian import Gaussian
from ase.io import read, write
from ase.geometry.analysis import Analysis
import glob


# add2executable
def extract_bonds(readable_file):
    """
    Guesses the bonds in a molecule according to
    ase.neighborlist.natural_cutoffs.

    Parameters
    ==========
    readable_file: str
        string to the configuration file in any of the ASE readable formats.

    Return
    ======
    (list) [#bonds x 2(int)] Bonds in the molecule.

    E.g.
    myutils extract_bonds optimization.xyz
    """
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
def diff_bonds(conf1, conf2, frozen_dofs='frozen_dofs.dat'):
    """
    This function returns the bonds that are in one conf but not in the other.

    Parameters
    ==========
    conf11: str
        first configuration file to be compared.
    conf2: str
        second configuration file to be compared.
    frozen_dofs: (optional)
        file with the frozen DOFs. If a rupture is obtained, it will be
        added to this file.

    Return
    ======
    (list) [#brocken_bonds x 2(int)] pair of brocken bonds

    Note
    ====
    The comparison is only in one direction, namely, this function does not
    find the bonds in conf2 that are not in conf1.
    """

    bonds1 = extract_bonds(conf1)
    bonds2 = extract_bonds(conf2)

    different_bonds = []
    for bond in bonds1:
        if bond not in bonds2:
            different_bonds.append(bond)

    with open(frozen_dofs, 'a') as fdofs:
        for bond in different_bonds:
            for index in bond:
                fdofs.write(str(index) + ' ')
            fdofs.write('F\n')
    return different_bonds


# add2executable
def conf2pdb(confile, pdbtemplate, pdboutput=None):
    """
    Transform a configuration in a file (xyz, log...) into a pdb using a pdb
    file as template.

    Parameters
    ==========
    confile: str
        path to the config file to be transformed to pdb.
    pdbtemplate: str
        path to the pdb template for the output.
    pdbfile: str (optional)
        name of trasnformed config file with pdb format. The default name is
        the same than the confile but with pdb extension.

    Return
    ======
    (str) The name of the output pdb file.

    Note
    ====
        the pdb file must contain the same atoms in the same order than the xyz
        file, this file only would change the coordinates.

    E.g.
    myutils optimization.log template.pdb
    myutils optimization.xyz template.pdb
    """
    if pdboutput is None:
        pdboutput = confile.split('.')[0] + '.pdb'

    atoms_ref = read(pdbtemplate)
    atoms_xyz = read(confile)
    assert len(atoms_ref) == len(atoms_xyz), "The number of atoms in the" + \
        " reference and the template does not coincide"
    assert \
        atoms_ref.get_chemical_symbols() == atoms_xyz.get_chemical_symbols(), \
        "The atoms in the reference and the template does not coincide"
    new_positions = atoms_xyz.positions.copy()

    atoms_ref.set_positions(new_positions)
    write(pdboutput, atoms_ref)

    return pdboutput


# add2executable
def all_xyz2pdb(template):
    """
    Transform all xyz files of the directory where it is executed into a pdb
    using a pdb file as template.

    Parameters
    ==========
    pdbtemplate: str
        path to the pdb template for the output.

    Return
    ======
    (str) The name of the outputs of each pdb file.

    Note
    ====
        the pdb file must contain the same atoms in the same order than the xyz
        file, this file only would change the coordinates.

    E.g.
    cd dir_with_xyz_files ; myutils all_xyz2pdb
    """
    configs = glob.glob('*.xyz')
    configs.sort()
    for config in configs:
        yield conf2pdb(config, template)


def all_hydrogen_atoms(mol):
    """"
    Finds the indexes of all the hydrogen atoms in the peptide from Atoms
    object.

    Parameters
    ==========
    mol: string to config file or ase.Atoms object
        ASE.Atoms object to extract the hydrogen indexes.

    Return
    ======
    (list) [#h_atoms(int)] Indexes corresponding to Hydrogen atoms in g09
        convention.
    """
    if type(mol) is str:
        mol = read(mol)
    indexes = np.where(mol.get_atomic_numbers() == 1)[0]

    return indexes + 1


# add2executable
def distance(file, index1, index2):
    """
    Computes the distance between two atoms in the last configuration
    of a trajectory file (e.g. .log file from gaussian).

    Parameters
    ==========
    file: str
        name of the file that contains the trajectory.
    index1: int
        index of the first atom to compute distances
    index2: int
        index of the second atom to compute distances

    Return
    ======
    (float)  Distance between atoms corresponding with atom with index1 and
    index2.

    E.g.
    myutils distance optimization.log 1 20
    """
    index1 = int(index1)
    index2 = int(index2)
    atoms = read(file)
    d = atoms.get_distance(index1, index2)

    return d


class MoleculeSetter:
    def __init__(self, atoms):
        self.atoms = atoms

    def rot_x(self, angle):
        """
        Rotation matrix around x axis.

        Parameters
        ==========
        angle: float[radians]
            angle to rotate around the x axis

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[1, 0, 0],
                      [0, c, -s],
                      [0, s, c]])
        return R

    def rot_y(self, angle):
        """
        Retuns the rotation matrix around y axis.

        Parameters
        ==========
        angle: float[radians]
            angle to rotate around the y axis

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, 0, s],
                      [0, 1, 0],
                      [-s, 0, c]])
        return R

    def rot_z(self, angle):
        """
        Retuns the rotation matrix around z axis.

        Parameters
        ==========
        angle: float[radians]
            angle to rotate around the z axis

        Return
        ======
        (numpy.array) [3float x 3float] Rotation matrix.
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, -s, 0],
                      [s, c, 0],
                      [0, 0, 1]])
        return R

    def align_axis(self, vector):
        """
        Apply the necessary rotations to set a vector aligned with positive x
        axis.

        Parameters
        ==========
        vector: array
            vector to be aligned.

        Return
        ======
        (numpy.array) [3float x 3float] Transformation matrix.
        """
        xyproj = vector.copy()
        xyproj[2] = 0
        phi = np.arcsin(vector[2]/np.linalg.norm(vector))
        theta = np.arccos(vector[0]/np.linalg.norm(xyproj))
        if vector[1] < 0:
            theta *= -1
        trans = np.dot(self.rot_y(phi), self.rot_z(-theta))
        return trans

    def align_plane(self, vector):
        """
        Rotation around x axis to set a vector in the xy plane.

        Parameters
        ==========
        vector: array
            vector to be rotated to be in the xy plane.

        Return
        ======
        (numpy.array) [3float x 3float] Transformation matrix.
        """
        reference = vector.copy()
        reference[0] = 0
        angle = np.arccos(reference[1]/np.linalg.norm(reference))
        if reference[2] < 0:
            angle *= -1
        return self.rot_x(-angle)

    def apply_trans(self, trans, indexes=None):
        """
        Apply a transformation to all vector positions of some atoms.

        Parameters
        ==========
        trans: array (3x3)
            transformation matrix to be applied to all atom positions.
        indexes: list or array (optional)
            indexes of the atoms to apply the transformation. Default None
            that means the transformation is applied to the positions of all
            the atoms.

        Return
        ======
        (numpy.array)[#natoms x 3float] new xyz positions of the N atoms. It
            changes the positions of the internal atoms object.
        """
        if indexes is None:
            indexes = list(range(len(self.atoms)))

        new_positions = []
        for i, atom in enumerate(self.atoms):
            if i in indexes:
                new_positions.append(np.dot(trans, atom.position))
            else:
                new_positions.append(atom.position)
        self.atoms.set_positions(new_positions)

        return new_positions

    def xy_alignment(self, index1, index2, index3=None, center=None):
        """
        Transform the positions of the atoms such that the atoms of indexes1
        and index2 are aligned in the x axis. The atom with index3 would be in
        the xy plane in case to be given.

        Parameters
        ==========
        index1 and index2: int
            indexes of the atoms to be aligned with the x-axis. The positive
            direction of x would go from atom 1 to atom 2.
        index 3: int (optional)
            The atom with index 3 would be in the xy plane in case to be given.
        Center: int (optional)
            It must be index1 or index2, that means the atom with this index
            will be placed in the origin. In case center=None (default), the
            origin would be in the geometrical center between atoms with index1
            and index2.

        Return
        ======
        (numpy.array)[#natoms x 3float] new xyz positions of the N atoms. It
            changes the positions of the internal atoms object.
        """
        # Move the origin
        if center == index1:
            center = self.atoms[center].position
        elif center == index2:
            center = self.atoms[center].position
            index2 = index1
            index1 = center
        else:
            center = (self.atoms[index1].position +
                      self.atoms[index2].position)/2

        self.atoms.set_positions(self.atoms.positions - center)
        # set index1 and index2 along x axis
        axis = self.atoms[index2].position
        self.apply_trans(self.align_axis(axis))
        if index3 is not None:
            third = self.atoms[index3].position
            self.apply_trans(self.align_plane(third))
        return self.atoms.positions

    def increase_distance(self, constraints, deltad):
        """
        increase the distance between two atoms by moving them and keeping the
        rest of the atoms in the same place.

        Parameters
        ==========
        constraints:
            constraints with the shape (n, 2), where n is the number of
            constraints and the first pair is the one to increase the
            distance.
        deltad: float
            amount to add to the distance between atoms.

        Return
        ======
        (ase.Atoms) Internal Atoms object with the corresponding modification.
        """
        d1norm = self.atoms.get_distance(constraints[0][0], constraints[0][1])
        self.atoms.set_distance(constraints[0][0], constraints[0][1],
                                d1norm+deltad)
        return self.atoms

    def increase_distance_with_constraints(self, constraints, deltad):
        """
        Take a configuration and increase the distance between two atoms by
        moving those atoms and all constraints containing them and keeping
        the rest of the atoms in the same place.

        Parameters
        ==========
        constraints:
            constraints with the shape (n, 2), where n is the number of
            constraints and the first pair is the one to increase the
            distance.
        deltad: float
            amount to add to the distance between atoms.

        Return
        ======
        (ase.Atoms) Internal Atoms object with the corresponding modification.
        """
        self.xy_alignment(constraints[0][0], constraints[0][1])
        left = [constraints[0][0]]
        right = [constraints[0][1]]
        len_right = 0
        len_left = 0
        while ((len_left != len(left)) and (len_right != len(right))):
            len_left = len(left)
            len_right = len(right)
            for cons in constraints[1:]:
                if ((cons[0] in left) and (cons[1] not in left)):
                    left.append(cons[1])
                elif ((cons[1] in left) and (cons[0] not in left)):
                    left.append(cons[0])
                elif ((cons[0] in right) and (cons[1] not in right)):
                    right.append(cons[1])
                elif ((cons[1] in right) and (cons[0] not in right)):
                    right.append(cons[0])
        print("These will be the atoms moved to the left: ", left)
        print("These will be the atoms moved to the right: ", right)

        new_positions = []
        for i, atom in enumerate(self.atoms):
            if i in right:
                new_positions.append(atom.position +
                                     np.array([deltad/2, 0, 0]))
            elif i in left:
                new_positions.append(atom.position +
                                     np.array([-deltad/2, 0, 0]))
            else:
                new_positions.append(atom.position)

        self.atoms.set_positions(new_positions)
        return self.atoms

    def scale_distance(self, constraints, deltad, index3=None):
        """
        Increase the distance between two atoms by aligning those atoms with
        the x axis and scaling the x-coordinate of all intermedia atoms.

        Parameters
        ==========
        constraints: list
            constraints with the shape (n, 2), where n is the number of
            constraints and the first pair is the one to increase the
            distance.
        deltad: float
            amount to add to the distance between atoms.

        Return
        ======
        (ase.Atoms) Internal Atoms object with the corresponding modification.
        """
        index1, index2 = constraints[0]
        d1norm = self.atoms.get_distance(index1, index2)
        # Move atom1 to the origin and rotate the molecule such that atom2 is
        # aligned with the +x axis:
        self.xy_alignment(index1, index2, center=index1)
        new_positions = [atom.position *
                         np.array([(d1norm + deltad)/d1norm, 1, 1])
                         for atom in self.atoms]
        self.atoms.set_positions(new_positions)
        if index3 is not None:
            third = self.atoms[index3].position
            self.apply_trans(self.align_plane(third))

        return new_positions

    def create_gaussian_input(self, out=None, charge=0, xc='bmk',
                              basis='6-31+g'):
        """
        Creates a g09 .com file without specifing the kind of calculus to run
        (optimization, ab-initio md, frequencies...). You would have to add it.

        Parameters
        ==========
        out: str (optional)
            name of the gaussian file (.com) without extension. Default
            chemical formula.
        charge: int (optional)
            charge of the molecule in electron units. Default 0.
        xc: str (optional)
            exchange correlation functional used in gaussian. Default bmk
        basis: str (optional)
            basis set used in gaussian. Default 6-31+g

        Return
        ======
        (ase.calculator.Gaussian) Calculator used to create the input.
        """
        if out is None:
            out = self.atoms.get_chemical_formula()

        calculator = Gaussian(label=out,
                              chk=out,
                              xc=xc,
                              basis=basis,
                              mult=1,
                              charge=int(charge))

        calculator.write_input(self.atoms)

        return calculator
