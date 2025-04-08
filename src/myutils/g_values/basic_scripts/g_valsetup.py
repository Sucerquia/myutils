from ase.io import read
from ase.geometry.analysis import Analysis
import numpy as np


def HFC_relevantA(atoms, radicals, depth=2):
    """
    Extract the H atoms of the neighbors up to certain depth (neighbors of
    neighbors). It also includes all O and N atoms.

    Parameters
    ==========
    atoms: ase.Atoms
        Atoms object containing the molecule.
    radicals: list
        posible positions of the radicals, implying that that neighborhood is
        important for the HFC.
    depth: int. Default=2
        number of times that it searches the neighbors of the neighbors.

    Return
    ======
    (np.array) Indexes of the H atoms belonging to the neighborhood of the
    radicals and the indexes of the N and O atoms. The indixes start with 1.
    """
    centers = radicals
    elements = np.array(atoms.get_chemical_symbols())
    ana = Analysis(atoms)
    relevant = []
    
    for j in range(depth):
        new_neighbors = []
        for i in centers:
            neighbors = np.array(ana.all_bonds[0][i])
            e_neighbors = elements[neighbors]
            hydrogens = e_neighbors == 'H'
            oxygens = e_neighbors == 'O'
            nitrogens = e_neighbors == 'N'
            condition = np.logical_or(np.logical_or(hydrogens,
                                                    oxygens),
                                      nitrogens)
            relevant.append(neighbors[condition])
            new_neighbors.append(neighbors[e_neighbors != 'H'])
        centers = [index for sublist in new_neighbors for index in sublist]
    relevant = np.array([index for sublist in relevant for index in sublist])
    return np.unique(relevant) + 1


# add2executable
def iHFC_fromxyz(file, radicals, depth):
    """
    Extract the H atoms of the neighbors up to certain depth (neighbors of
    neighbors). It also includes all O and N atoms. Those atoms are the
    relevant ones for HyperFine Calculations (HFC). This function receives
    strings as inputs.

    Parameters
    ==========
    file: str
        path to the molecule that you want to analyse. ASE must be able to read
        this file.
    radicals: str
        index of the atoms of posible location of the radicals, implying that
        that neighborhoods of this atoms are important for the HFC. It starts
        from 0 and it should
        have the shape of a list, f.e. '[0, 7]'.
    depth: str
        number of times that it searches the neighbors of the neighbors.

    Return
    ======
    (np.array) Indexes of the H atoms belonging to the neighborhood of the
    radicals and the indexes of the N and O atoms. It starts from 1.
    """
    atoms = read(file)
    radicals = eval(radicals)
    depth = int(depth)

    return HFC_relevantA(atoms, radicals, depth)


# add2executable
def orca_opt_inp(xc, base, processors, multiplicity, charge, xyzFile):
    if xc == '':
        xc = 'B3LYP'
    if base == '':
        base = 'EPR-II'
    if processors == '':
        processors = '16'
    if multiplicity == '':
        multiplicity = '1'
    if charge == '':
        charge = '0'
    if xyzFile == '':
        xyzFile = 'model.xyz'

    print(f"! {xc} {base} OPT")
    print(f"%pal nprocs {processors} end")
    print(f"*XYZFile {charge} {multiplicity} {xyzFile}")


# add2executable
def orca_epr_inp(xc, base, processors, multiplicity, charge, xyzFile, nuclei):
    if xc == '':
        xc = 'B3LYP'
    if base == '':
        base = 'EPR-II'
    if processors == '':
        processors = '16'
    if multiplicity == '':
        multiplicity = '2'
    if charge == '':
        charge = '0'
    if xyzFile == '':
        xyzFile = 'opt_EPRII.xyz'

    print(f"! {xc} {base} AUTOAUX")
    print(f"%pal nprocs {processors} end")
    print(f"*XYZFile {charge} {multiplicity} {xyzFile}")
    print("%EPRNMR\n        GTENSOR   TRUE")
    if nuclei != '':
        all_nuclei = nuclei.split('-')
        for nuclei in all_nuclei:
            print("        NUCLEI    = " + nuclei + " {SHIFT, AISO, ADIP, AORB}")
    print("        ORI       GIAO\nEND")