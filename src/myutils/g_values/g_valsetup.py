# %%
from ase.io import read
from ase.geometry.analysis import Analysis
import numpy as np

import networkx as nx


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


def get_connectivity(atoms, engine):
    if engine == 'vmd':
        from ase.io import write
        from myutils.ase_utils.molecules import vmd_connectivity
        import os


        write('tmp_Mview.xyz', atoms)
        connectivity = vmd_connectivity('tmp_Mview.xyz')
        os.remove('tmp_Mview.xyz')

    elif engine == 'ase':
        from ase.neighborlist import natural_cutoffs, NeighborList


        cutoffs = natural_cutoffs(atoms)
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)
        connectivity = [nl.get_neighbors(i)[0] for i in range(len(atoms))]
    
    return connectivity


def ase_to_nx_graph(atoms, engine='vmd'):
    connectivity  = get_connectivity(atoms, engine)

    G = nx.Graph()
    for i, atom in enumerate(atoms):
        G.add_node(i, atomic_number=atom.number)

    for i in range(len(atoms)):
        indices = connectivity[i]
        for j in indices:
            if not G.has_edge(i, j):
                G.add_edge(i, j)

    return G

def get_isomorphic_permutation(nx_graph1, nx_graph2):
    """
    Returns a permutation list to reorder the atoms of `graph2` to match `graph1` 
    based on their isomorphic structure such that e.g.:

    permutation = get_isomorphic_permutation(g1, g2)
    g2.atomic_number[:,permutation] = g1.atomic_number

    Parameters:
    - graph1: A DGLGraph object representing the first molecular graph.
    - graph2: A DGLGraph object representing the second molecular graph.

    Returns:
    - List[int]: A list of indices representing the permutation required to align graph1 to graph2.
    """

    # Convert DGLGraphs to NetworkX graphs for isomorphism checking
    def node_match(n1, n2):
        return n1['atomic_number'] == n2['atomic_number']
        
    # Find the isomorphism mapping between graph1 and graph2 nodes
    gm = nx.isomorphism.GraphMatcher(nx_graph1, nx_graph2, node_match=node_match)

    if gm.is_isomorphic():
        # Extract the node correspondence mapping from graph1 to graph2
        mapping = gm.mapping
        # Generate the permutation list based on the mapping
        permutation = [mapping[i] for i in range(len(nx_graph1))]
        return permutation
    else:
        raise ValueError("Graphs are not isomorphic")

def get_permutation_of_heavy_atoms(heavy_ref_mol, heavy_mch_mol, engine='vmd'):
    graph_ref = ase_to_nx_graph(heavy_ref_mol, engine=engine)
    graph_2ma = ase_to_nx_graph(heavy_mch_mol, engine=engine)
    perm = get_isomorphic_permutation(graph_ref, graph_2ma)

    return perm

# add2executable
def rad_loc(ref_mol_xyz: str, rad_mol_xyz: str, engine='vmd') -> int:
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

    ref_mol = read(ref_mol_xyz)
    con_ref = get_connectivity(ref_mol, engine)
    elements = np.array(ref_mol.get_chemical_symbols()).copy()
    non_H_ref = np.where(elements != 'H')[0]
    heavy_ref_mol = ref_mol[elements != 'H']
    
    mch_mol = read(rad_mol_xyz)
    con_mch = get_connectivity(mch_mol, engine)
    elements = np.array(mch_mol.get_chemical_symbols()).copy()
    non_H_mch = np.where(elements != 'H')[0]
    heavy_mch_mol = mch_mol[elements != 'H']

    perm = get_permutation_of_heavy_atoms(heavy_ref_mol,
                                          heavy_mch_mol,
                                          engine=engine)
    
    for i, j in zip(non_H_ref, non_H_mch[perm]):
        if len(con_ref[i]) != len(con_mch[j]):
            return j + 1


# add2executable
def heavya_idx_from_ref(ref_mol_xyz: str, rad_mol_xyz: str,
                        index: int, engine: str='vmd') -> np.ndarray:
    """
    Find the index of the atom where a hydrogen atom was attached.

    Parametes
    =========
    ref_mol: str
        path to the file describing the molecule with the hydrogen included.
    redical: str
        path to the file describing the molecule without the hydrogen.
    index: int
        index of the heavy atom in the reference molecule to match in the new
        molecule. Starting from 1.

    Return
    ======
    (int) index of the atom where the hydrogen was attached. Index counting
    starting from 1.
    """

    ref_mol = read(ref_mol_xyz)
    elements = np.array(ref_mol.get_chemical_symbols()).copy()
    non_H_ref = np.where(elements != 'H')[0]
    heavy_ref_mol = ref_mol[elements != 'H']
    
    mch_mol = read(rad_mol_xyz)
    elements = np.array(mch_mol.get_chemical_symbols()).copy()
    non_H_mch = np.where(elements != 'H')[0]
    heavy_mch_mol = mch_mol[elements != 'H']

    perm = get_permutation_of_heavy_atoms(heavy_ref_mol, heavy_mch_mol,
                                          engine=engine)
    
    idx_in_reduced_ref = np.where(non_H_ref == int(index) - 1)[0]
    if len(idx_in_reduced_ref) != 1:
        raise ValueError(f"{index} seems not to be a heavy atom or it is not" +
                         " part of the molecule to match.")

    return non_H_mch[perm][idx_in_reduced_ref][0] + 1


def rad_loc_deprecated(ref_mol: str, radical: str) -> np.ndarray:
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
