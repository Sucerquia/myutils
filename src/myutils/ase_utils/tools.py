import numpy as np
from ase.calculators.gaussian import Gaussian


class MoleculeSetter:
    def __init__(self, atoms):
        self.atoms = atoms

    def rot_x(self, angle):
        """
        Retuns the rotation matrix around x axis.

        Parameters
        ==========
        angle: float[radians]
            angle to rotate around the x axis
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
        """
        c = np.cos(angle)
        s = np.sin(angle)
        R = np.array([[c, -s, 0],
                      [s, c, 0],
                      [0, 0, 1]])
        return R

    def align_axis(self, vector):
        """
        Apply the necessary rotations to set a
        vector aligned with positive x axis.

        Parameters
        ==========
        vector: array
            vector tp be aligned.

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
        """
        reference = vector.copy()
        reference[0] = 0
        angle = np.arccos(reference[1]/np.linalg.norm(reference))
        if reference[2] < 0:
            angle *= -1
        return self.rot_x(-angle)

    def apply_trans(self, trans):
        """
        Apply a transformation to all vector positions of the atoms object.

        Parameters
        ==========
        trans: array (3x3)
            transformation matrix to be applied to all atoms positions.
        """
        new_positions = [np.dot(trans, atom.position) for atom in self.atoms]
        self.atoms.set_positions(new_positions)
        return new_positions

    def xy_alignment(self, index1, index2, index3=None, center=None):
        """
        Transforme the positions of the atoms such that the atoms of indexes 1
        and 2 are aligned in the x axis. The atom with index 3 would be in the
        xy plane in case to be given.

        Parameters
        ==========
        index 1 and 2: int
            indexes of the atoms to be aligned with the x-axis. The positive
            direction of x would go from atom 1 to atom 2.

        index 3: int (optional)
            The atom with index 3 would be in the xy plane in case to be given.

        Center: int (optional)
            It must be index1 or index2, that means the atom with this index
            will be placed in the origin. In case center=None (default), the
            origin would be in the geometrical center between atoms with index1
            and index2.
        """
        # Move the origin
        if center is None:
            center = (self.atoms[index1].position +
                      self.atoms[index2].position)/2
        else:
            center = self.atoms[center].position
        self.atoms.set_positions(self.atoms.positions - center)
        # set index1 and index2 along x axis
        axis = self.atoms[index2].position
        self.apply_trans(self.align_axis(axis))
        if index3 is not None:
            third = self.atoms[index3].position
            self.apply_trans(self.align_plane(third))
        return self.atoms.positions

    def increase_distance(self, index1, index2, deltad):
        """
        Increase the distance between two atoms keeping
        the same midpoint.

        Parameters
        ==========
        index1: int
            index of one of the atoms to increase the distance.
        index2: int
            index of the other atom to increase the distance.
        deltad: float
            amount to add to the distance between atoms.
        """
        d1norm = self.atoms.get_distance(index1, index2)
        self.atoms.set_distance(index1, index2, d1norm+deltad)
        return self.atoms

    def increase_distance_with_constrain(self, index1, index2, deltad,
                                         constrains):
        """
        Increase the distance between two atoms keeping
        the same midpoint.

        Parameters
        ==========
        index1: int
            index of one of the atoms to increase the distance, ASE notation.
        index2: int
            index of the other atom to increase the distance, ASE notation.
        deltad: float
            amount to add to the distance between atoms.
        constrains:
            pairs of atoms to keep the same distance, ASE notation. This
            parameter has to be defined as:
            [[cons1_ind1, cons1_ind2], [cons2_ind1, cons2_ind2], ... ]
            in g09 index notation.
        """
        left = [index1]
        right = [index2]
        len_right = 0
        len_left = 0
        while ((len_left != len(left)) and (len_right != len(right))):
            len_left = len(left)
            len_right = len(right)
            for cons in constrains:
                if ((cons[0] in left) and (cons[1] not in left)):
                    left.append(cons[1])
                elif ((cons[1] in left) and (cons[0] not in left)):
                    left.append(cons[0])
                elif ((cons[0] in right) and (cons[1] not in right)):
                    right.append(cons[1])
                elif ((cons[1] in right) and (cons[0] not in right)):
                    right.append(cons[0])
        print("this will be the atoms moved to the left: ", left)
        print("this will be the atoms moved to the right: ", right)

        new_positions =[]
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

    def scale_distance(self, index1, index2, deltad, index3=None):
        """
        Increase the distance between two atoms by aligning those atoms with
        the x axis and scaling the x-coordinate of all intermedia atoms.

        Parameters
        ==========
        index1: int
            index of one of the atoms to increase the distance.
        index2: int
            index of the other atom to increase the distance.
        deltad: float
            amount to add to the distance between atoms.
        """
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

    # Not working yet
    def scale_w_cons(self, index1, index2, deltad, constrains):
        """
        Increase the distance between two atoms by aligning those atoms with
        the x axis and scaling the x-coordinate of all intermedia atoms.

        Parameters
        ==========
        index1: int
            index of one of the atoms to increase the distance.
        index2: int
            index of the other atom to increase the distance.
        deltad: float
            amount to add to the distance between atoms.
        constrains:
            pairs of atoms to keep the same distance. with the shape:
            [[cons1_ind1, cons1_ind2], [cons2_ind1, cons2_ind2], ... ]
            in g09 index notation.
        """
        d1norm = self.atoms.get_distance(index1, index2)
        # Move atom1 to the origin and rotate the molecule such that atom2 is
        # aligned with the +x axis:
        self.xy_alignment(index1, index2, center=index1)

        frozen = []
        for i in constrains.flatten():
            if i not in frozen:
                frozen.append(i)

        new_positions=[]
        for i, atom in enumerate(self.atoms):
            if i + 1 not in frozen:
                new_positions.append(atom.position *
                                 np.array([(d1norm + deltad)/d1norm, 1, 1]))
        self.atoms.set_positions(new_positions)

        return new_positions

    def create_gaussian_input(self, out=None, charge=0, xc='bmk',
                              basis='6-31+g'):
        """
        Creates a g09 .com file without specifing the kind of calculus to run
        (optimization, ab-initio md, frequencies...).

        Parameters
        ==========
        out: str
            name of the gaussian file (.com) without extension.
        charge: int [e]
            charge of the molecule in electron units. Default 0.
        xc: str
            exchange correlation functional used in gaussian. Default bmk
        basis: str
            basis set used in gaussian. Default 6-31+g
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
