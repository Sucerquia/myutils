import numpy as np
from myutils.peptides import PepSetter
from myutils.miscellaneous import output_terminal
from ase.io import read
import sys
from pytest import approx


def extract_dofs(indexes, atoms):
    distance = atoms.get_distance(*indexes[:2])
    angle = atoms.get_angle(*indexes[:3])
    dihedral = atoms.get_dihedral(*indexes)
    while dihedral > 180:
        dihedral -= 360

    return distance, angle, dihedral


def test_dofs(atoms, indexes, file):
    dist, angl, dihe = extract_dofs(np.array(indexes) - 1, atoms)
    old_dist = float(output_terminal(f'sed -n "/R{indexes[0]}=/p" {file}' +
                                      ' | cut -d = -f 2', print_output=False))
    old_angl = float(output_terminal(f'sed -n "/A{indexes[0]}=/p" {file}' +
                                      ' | cut -d = -f 2', print_output=False))
    old_dihe = float(output_terminal(f'sed -n "/D{indexes[0]}=/p" {file}' +
                                      ' | cut -d = -f 2', print_output=False))

    assert old_dist == approx(dist, abs=3e-2), f"R{indexes[0]}({old_dist}) does not " + \
        f"correspond to the expected from the xyz file({dist})"
    assert old_angl == approx(angl, abs=5e-1), f"A{indexes[0]}({old_angl}) does not " + \
        f"correspond to the expected from the xyz file({angl})"
    assert old_dihe == approx(dihe, abs=5e-1), f"D{indexes[0]}({old_dihe}) does not " + \
        f"correspond to the expected from the xyz file({dihe})"


def def_line(indexes, element):
    a1, a2, a3, a4 = indexes
    return f'{element},{a2},R{a1},{a3},A{a1},{a4},D{a1}'


def change_def(old_i, new_i, element, atoms, file):
    # Create new and old lines
    old_line = def_line(old_i, element)
    new_line = def_line(new_i, element)

    # Find variables
    new_i = np.array(new_i) - 1
    old_i = np.array(old_i) - 1
    dist, angl, dihe = extract_dofs(new_i, atoms)
    
    # Change value
    output_terminal(f'sed -i "s/{old_line}/{new_line}/g" {file}')
    output_terminal(f'sed -i "/R{old_i[0]}=/c\ R{old_i[0]}={dist}" {file}')
    output_terminal(f'sed -i "/A{old_i[0]}=/c\ A{old_i[0]}={angl}" {file}')
    output_terminal(f'sed -i "/D{old_i[0]}=/c\ D{old_i[0]}={dihe}" {file}')


# add2executable
def change_prolines_dofs(comfile, molecule, pdb_template):
    pep_set = PepSetter(pdb_template)
    atoms = read(molecule)
    pros_i = np.where(np.array(list(pep_set.amino_name.values())) == 'PRO ')[0] + 1

    for i_pro in pros_i:
        amino = pep_set.amino_info[i_pro]
        Ca_i = amino['CA']
        Cb_i = amino['CB']
        Cg_i = amino['CG']
        Cd_i = amino['CD']
        HG1_i = amino['2HG']
        HG2_i = amino['3HG']
        N_i = amino['N']

        # test that xyz does does correspond to the com file
        test_dofs(atoms,
                  np.array([Cg_i, Cb_i, Ca_i, N_i]),
                  comfile)
        test_dofs(atoms,
                  np.array([HG1_i, Cg_i, Cb_i, Ca_i]),
                  comfile)
        test_dofs(atoms,
                  np.array([HG2_i, Cg_i, Cb_i, HG1_i]),
                  comfile)

        # change Cg dofs
        change_def(np.array([Cg_i, Cb_i, Ca_i, N_i]),
                   np.array([Cg_i, Cd_i, N_i, Ca_i]),
                   'C', atoms, comfile)

        # change Gg1 dofs
        change_def(np.array([HG1_i, Cg_i, Cb_i, Ca_i]),
                   np.array([HG1_i, Cg_i, Cd_i, N_i]),
                   'H', atoms, comfile)

        # change Gg1 dofs
        change_def(np.array([HG2_i, Cg_i, Cb_i, HG1_i]),
                   np.array([HG2_i, Cg_i, Cd_i, N_i]),
                   'H', atoms, comfile)
        
        output_terminal(f'myutils switch_atoms_in_com {Cg_i} {Cd_i} {comfile}')
