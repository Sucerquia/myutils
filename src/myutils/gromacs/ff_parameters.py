from openmm.app import PDBFile, ForceField
from grappa.utils.model_loading_utils import model_from_tag
from grappa.data import Parameters
from grappa.data import Molecule


def create_bondpars_file(outfile, bond_idxs, bond_eqs, bond_ks):
    """
    Creates a file containing the bond parameters.

    Parameters
    ==========
    outfile: str
        path to the outputfile. Tipically a dat file.

    Return
    ======
    (tuple) bond_idxs, bond_eqs, bond_ks
    """
    with open(outfile, "w") as grapFile:
        grapFile.write("# index1 index2 k_val[kcal/mol/A^2] r_0[A]")
        for i in range(len(bond_eqs)):
            grapFile.write(str(bond_idxs[i][0]) + "\t"
                        + str(bond_idxs[i][1]) + "\t"
                        + str(bond_ks[i]) + "\t"
                        + str(bond_eqs[i]) + "\n")

    return bond_idxs, bond_eqs, bond_ks


def grappa_bond_pars(pdb):
    openmm_topology = PDBFile(pdb).getTopology()
    openmm_system = ForceField(
                               'amber99sbildn.xml'
                               ).createSystem(openmm_topology)

    molecule = Molecule.from_openmm_system(openmm_system= openmm_system,
                                           openmm_topology=openmm_topology)

    # molecule -> dgl(g) -> parameters
    g = molecule.to_dgl()
    model = model_from_tag('grappa-1.3')
    g = model(g)

    # zero-based idxs of the atoms in the respective bond (corresponding to the
    # order in the openmm topology, i.e. the order in the PDB file)
    # shape (n_bonds, 2)
    bond_idxs = g.nodes['n2'].data['idxs'].detach().cpu().numpy()

    # bond lengths in angstrom
    # shape (n_bonds,)
    bond_eqs = g.nodes['n2'].data['eq'].detach().cpu().numpy()

    # bond force constants in kcal/mol/angstrom^2
    # shape (n_bonds,)
    bond_ks = g.nodes['n2'].data['k'].detach().cpu().numpy()

    return bond_idxs, bond_eqs, bond_ks


# add2executable
def create_grappa_data(pdb, outfile):
    """
    Create the file containing the parameters according to grappa.

    Parameters
    ==========
    pdb: str
        path to pdb file containing the basic structure. (it is assumed that
        the topollogy remain unchanged)
    outfile: str
        path to the outputfile. Tipically a dat file.

    Return
    ======
    (tuple) bond_idxs, bond_eqs, bond_ks  created by grappa.
    """
    bond_idxs, bond_eqs, bond_ks = grappa_bond_pars(pdb)
    create_bondpars_file(outfile, bond_idxs, bond_eqs, bond_ks)
    return bond_idxs, bond_eqs, bond_ks


def amber_bond_pars(pdb):
    openmm_topology = PDBFile(pdb).getTopology()
    openmm_system = ForceField(
                               'amber99sbildn.xml'
                               ).createSystem(openmm_topology)

    molecule = Molecule.from_openmm_system(openmm_system= openmm_system,
                                           openmm_topology=openmm_topology)
    
    params = Parameters.from_openmm_system(openmm_system=openmm_system, mol=molecule)

    g = molecule.to_dgl()
    g = params.write_to_dgl(g, suffix='_ref')

    # repeat the procedure from above:
    bond_eqs = g.nodes['n2'].data['eq_ref'].detach().cpu().numpy()
    bond_ks = g.nodes['n2'].data['k_ref'].detach().cpu().numpy()
    bond_idxs = g.nodes['n2'].data['idxs'].detach().cpu().numpy()

    return bond_idxs, bond_eqs, bond_ks


# add2executable
def create_amber_data(pdb, outfile):
    """
    Create the file containing the parameters according to grappa.

    Parameters
    ==========
    pdb: str
        path to pdb file containing the basic structure. (it is assumed that
        the topollogy remain unchanged)
    outfile: str
        path to the outputfile. Tipically a dat file.

    Return
    ======
    (tuple) bond_idxs, bond_eqs, bond_ks created by amber.
    """
    bond_idxs, bond_eqs, bond_ks = amber_bond_pars(pdb)

    create_bondpars_file(outfile, bond_idxs, bond_eqs, bond_ks)

    return bond_idxs, bond_eqs, bond_ks