from myutils.ase_utils.molecules import MoleculeSetter
from ase.geometry.analysis import Analysis
from ase.io import read, write
import numpy as np
import glob


# add2executable
def change_distance(inp, out, file_cons, deltad, charge, method):
    """
    Take a configuration and increase the distance between two atoms. With the
    new structure, it creates a gaussian file without specifing the kind of
    calculus to run (optimization, ab-initio md, frequencies...).

    Parameters
    ==========
    inp: str
        input file containing the structure to modify.
    out: str
        name of the gaussian file (.com) without extension.
    file_cons: str
        file with the constraints sorted in the first two columns. The first
        pair is the one to change the distance.
    deltad: float
        amount to add to the distance between atoms.
    charge: int
        charge in e to create the g09 input file. The multiplicity is assumed.
        to be one.
    method: str
        method defined in myutils.ase_utils.tools. So far, the methods already
        implemented are: 'scale_distance', 'increase_distance',
        'increase_distance_with_constraints'.
    """
    methods = ['scale_distance', 'increase_distance',
               'increase_distance_with_constraints']
    if method not in methods:
        raise ValueError("Non-recognized stretching method. To see the " +
                         "options, check 'myutils change_distance -h'")
    # -1 to transform into python convention
    cons = np.loadtxt(file_cons, usecols=[0, 1], dtype=int) - 1
    if type(cons[0]) is np.int64:
        cons = np.array([list(cons)])
    deltad = float(deltad)

    # Read previus file
    atoms = read(inp)
    manipulator = MoleculeSetter(atoms)
    manipulator.xy_alignment(cons[0][0], cons[0][1], center=cons[0][0])
    eval(f'manipulator.{method}(cons, deltad)')
    manipulator.create_gaussian_input(out=out, charge=charge)

    return f"{out}.com file created"


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
def all_xyz2pdb(template, output_patern=None, xyzdir=''):
    """
    Transform all xyz files of the directory where it is executed into a pdb
    using a pdb file as template.

    Parameters
    ==========
    pdbtemplate: str
        path to the pdb template for the output.
    output_patern: str (optional)
        the name of the output will be <this string>-<n>.pdb, where is is an
        increasing index, from 1 to the number of xyz files.

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
    configs = glob.glob(xyzdir + '*.xyz')
    configs.sort()
    n=1
    outfiles = []
    for config in configs:
        if output_patern is None:
            pdboutput = config.split('.')[0] + '.pdb'
        else:
            pdboutput = output_patern + f'-{n}.pdb'
            print(pdboutput)
            n += 1
        outfiles.append(conf2pdb(config, template, pdboutput=pdboutput))
    return outfiles


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
