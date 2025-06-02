from ase.io import read, write
from ase.geometry.analysis import Analysis
from ase import Atoms
import numpy as np


# add2executable
def rad_loc(ref_mol: str, radical: str) -> np.ndarray:
    """
    Find the index of the atom where a hydrogen atom was attached.

    Parametes
    =========
    ref_mol: str
        path to the file describing the molecule with the hydrogen included.
    redical: str
        path to the file describing the molecule without the hydrogen.

    Return
    ======
    (int) index of the atom where the hydrogen was attached. Index counting
    starting from 1.
    """
    withrad = read(radical)
    ana_rad = Analysis(withrad)
    cone_rad = ana_rad.all_bonds
    reference = read(ref_mol)
    Hs = np.where(np.array(reference.get_chemical_symbols()) == 'H')[0]

    for h in Hs:
        reference = read(ref_mol)
        del reference[h]

        ana = Analysis(reference)
        cone = ana.all_bonds

        if cone == cone_rad:
            reference = read(ref_mol)
            ana = Analysis(reference)
            cone = ana.all_bonds
            return cone[0][h][0] + 1
    raise ValueError("Conection to H not found")


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
        important for the HFC. Indexes starting from 1.
    depth: int. Default=2
        number of times that it searches the neighbors of the neighbors.

    Return
    ======
    (np.array) Indexes of the H atoms belonging to the neighborhood of the
    radicals and the indexes of the N and O atoms. The indices start with 1.
    """
    centers = np.array(radicals) - 1
    elements = np.array(atoms.get_chemical_symbols())
    ana = Analysis(atoms)
    relevant = []
    
    for _ in range(depth):
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
    
    # Add last hydrogens
    centers = np.array([index for sublist in relevant for index in sublist])
    for i in centers:
        neighbors = np.array(ana.all_bonds[0][i])
        e_neighbors = elements[neighbors]
        hydrogens = e_neighbors == 'H'
        relevant.append(neighbors[hydrogens])
    relevant = np.array([index for sublist in relevant for index in sublist])
    relevant = np.unique(relevant)
    separated_relevant = {}
    for i, element  in enumerate(elements[relevant]):
        if element not in list(separated_relevant.keys()):
            separated_relevant[str(element)] = []
        separated_relevant[element].append(int(relevant[i] + 1))

    return separated_relevant


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
def ext_xyz_from_npz(npz_file):
    data = np.load(npz_file)

    for i, xyz in enumerate(data['xyz']):
        atoms = Atoms(numbers=data['atomic_numbers'][i],
                      positions=xyz)
        name = npz_file[:-4] + f'_{i:03}' + '.xyz'

        comment = f'total_charge: {data["total_charge"]}; ' + \
                   f'multiplicity: 2; ' + \
                   f'heavy_atom_missing_idxs: ' + \
                   f'{data["heavy_atom_missing_idxs"][i]}'
        if 'atom_chargerelevant_idx' in data and len(data["atom_chargerelevant_idx"][i]) != 0:
            charged = ', '.join([str(j) for j in
                                 data["atom_chargerelevant_idx"][i]])
            comment += '; atom_chargerelevant_idx: ' + \
                       f'{charged}'
        if 'source_names' in data:
            comment += f'; source_name: {data["source_names"][i]}'
        write(name, atoms, comment=comment)
    
    info = {key: data[key] for key in data.files if key != 'original_xyz'}
    data.close()
    return info
   
