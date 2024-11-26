import numpy as np
from ase.io import read
import sys
import matplotlib.pyplot as plt
from myutils.miscellaneous import output_terminal
from ase.visualize import view
from myutils.ase_utils.molecules import MoleculeSetter, Alignment
from myutils.peptides import PepSetter
from ase import Atoms
from ase.io import write
import glob


# add2executable
def info_from_opt(pep):
    """
    Extract configurations
    Parameters
    ==========
    """
    # read configurations
    atoms = read(pep + "-optext.log", index=':')

    # read energies
    energies = output_terminal("grep 'SCF Done:' " + pep +
                               "-optext.log | awk '{print $5}'",
                               print_output=False)
    energies = np.array(energies.split('\n')[:-1], dtype=float)

    # read amino acids info, recognize extreme indexes
    ps = PepSetter(f'{pep}-stretched00.pdb')
    ind1 = ps.amino_info[1]['CH3'] - 1
    ind2 = ps.amino_info[5]['CH3'] - 1

    # choose third atom for orientation.
    if ps.amino_name[3] != 'GLY':  # Glycine does not have CB
        ind3 = ps.amino_info[3]['CB'] - 1
    elif ps.amino_name[2] != 'GLY':
        ind3 = ps.amino_info[2]['CB'] - 1
    elif ps.amino_name[4] != 'GLY':
        ind3 = ps.amino_info[4]['CB'] - 1
    else: # In case of 3 Glycine
        ind3 = ps.amino_info[3]['CA'] - 1

    distances = [conf.get_distance(ind1-1, ind2-1) for conf in atoms]
    energies2 = energies

    # remove configurations that goes up in energy. keep those that goes
    # down only. assuming local optimization
    i=0
    while True:
        de = energies[1:] - energies[:-1]
        toremove = np.where(de > 0)[0]
        if len(toremove) == 0:
            break

        for index in toremove[::-1]:
            atoms.pop(index + 1)
        energies = np.delete(energies, toremove+1)
        i += 1

    # align all the structures in the same plane. This guarantee that
    # the average of two structures is the intermedia structure between them
    atoms  = [Alignment.align_with_components(conf) for conf in atoms]
    for conf in atoms:
        ms = MoleculeSetter(conf)
        ms.xy_alignment(ind1, ind2, ind3)

    # TODO: change from here to consider continuos DOFs instead of continuos
    # distance. This function can be unified in this and the next one.

    # now make the trajectory continuos. if the distance between extremes is
    # larger than 0.2A, configurations with the intermedia distances are
    # created.
    new_set = [atoms[0]]
    i = 1
    while i < len(atoms):
        conf = atoms[i]
        di = new_set[-1].get_distance(ind1, ind2)
        df = conf.get_distance(ind1, ind2)
        deltad = df - di
        if abs(deltad) > 0.2:
            n_intermedia = int(abs(deltad / 0.2))
            fraction = 1/(n_intermedia + 1)
            inbetween = []
            for n in range(1, n_intermedia + 1):
                at = Atoms(conf.get_chemical_symbols(),
                           positions = (1 - fraction * n) *
                                       new_set[-1].positions +
                                       (fraction * n) * conf.positions)
                inbetween.append(at)
            new_set.extend(inbetween)
        new_set.append(conf)
        i += 1

    # in case of last configurations does not belong to new_set, it is added
    if np.any(atoms[-1].positions != new_set[-1].positions):
        new_set.append(atoms[-1])

    # write all the trajectory
    for i, atoms in enumerate(new_set[::-1]):
        write('{}-forces{:03d}.xyz'.format(pep, i), atoms)

# add2executable
def reduce_structs_pre(dir):
    """Check all the *-dofs.dat files and remove those files that represent
    unnrelevant changes"""
    all_files = glob.glob(f"{dir}/*-dofs.dat")
    all_files.sort()

    # Next variable is used later to guarantee same dofs and as DOFs definition
    dofs_ref = np.loadtxt(all_files[0],
                          delimiter='=',
                          comments='      Variables:',
                          usecols=0,
                          dtype=str)

    # Number of distances in the DOFs
    nrs = len([r for r in dofs_ref if r[1] == 'R'])

    all_dofs = []
    for file in all_files:
        # check that number of each kind of dofs is the same
        dofs = np.loadtxt(file,
                          delimiter='=',
                          comments='      Variables:',
                          usecols=0,
                          dtype=str)
        assert (dofs == dofs_ref).all(), \
            f"{file} has different dofs than {all_files[0]}"
        dofs = np.loadtxt(file,
                          delimiter='=',
                          comments='      Variables:',
                          usecols=1)
        all_dofs.append(dofs)
    all_dofs = np.array(all_dofs)

    lowest_dist = 0.01
    lowest_angl = 5
    i = 0
    j = 1
    d_ij = (all_dofs[j] - all_dofs[j + 1:])
    d_ij[:, :nrs] *= 10

    # make angles berween -180 to 180
    condition = d_ij[:, nrs:] < -180
    while condition.any():
        d_ij[:, nrs:][condition] += 360
        condition = d_ij[:, nrs:] < -180
    condition = d_ij[:, nrs:] > 180
    while condition.any():
        d_ij[:, nrs:][condition] -= 360
        condition = d_ij[:, nrs:] > 180

    # reduces by 2 the order of the distances
    d_ij[:, nrs:] *= 1e-2
    maxis_d = np.amax(d_ij[:, :nrs], axis=1)
    maxis_a = np.amax(d_ij[:, nrs:], axis=1)
    with open("file.dat", "w") as f:
        for i, j in zip(maxis_d, maxis_a):
            f.write(f"{i}\t {j}\n")
    """
        maxis = np.amax(d_ij, axis=0)
        close = np.isclose(d_ij, 0, atol=1e-3)
        print("close: ", close.shape)
        by_struc = np.all(close, axis=1)
        print("bs: ", by_struc.shape)
        same = np.where(by_struc)[0]
        print(same.shape)
        if len(same) != 0:
            print(j, np.where(by_struc))
    """
    new_set = []

    # first jump to repeated structures
    j = 0  # structures
    while j < len(all_files) - 1:
        new_set.append(j)
        d_ij = (all_dofs[j] - all_dofs[j + 1:])
        # rescale distances and angles to have them in the same order of
        # magnitude. I could also evaluate the approximation for distances and then
        # for angles and use logicaland.
        condition = d_ij[:, nrs:] < -180
        while condition.any():
            d_ij[:, nrs:][condition] += 360
            condition = d_ij[:, nrs:] < -180
        condition = d_ij[:, nrs:] > 180
        while condition.any():
            d_ij[:, nrs:][condition] -= 360
            condition = d_ij[:, nrs:] > 180
        d_ij[:, nrs:] *= 1e-3
        d_ij = abs(d_ij)
        close = np.isclose(d_ij, 0, atol=1e-3)
        by_struc = np.all(close, axis=1)
        same = np.where(by_struc)[0] + j
        if len(same) != 0:
            j = same[-1]
        j += 1
    print("bystr: ", new_set)

    new_dofs = all_dofs[new_set]
    deltas = new_dofs[1:] - new_dofs[:-1]
    maxis_d = np.amax(deltas[:, :nrs], axis=1)
    maxis_a = np.amax(deltas[:, nrs:], axis=1)
    condition = maxis_a < -180
    while condition.any():
        maxis_a[condition] += 360
        condition = maxis_a < -180
    condition = maxis_a > 180
    while condition.any():
        maxis_a[condition] -= 360
        condition = maxis_a > 180
    with open("file.dat", "w") as f:
        for i, j in zip(maxis_d, maxis_a):
            f.write(f"{i}\t {j}\n")
    with open("file2.dat", "w") as f:
        for j in new_set:
            f.write(f"{j}\n")


# add2executable
def reduce_structs(dir):
    """
    Check all the *-dofs.dat files and remove those files that represent
    irrelevant changes. It does not creates intermedias.
    """
    # find dofs files
    all_files = glob.glob(f"{dir}/*-dofs.dat")
    all_files.sort()

    # extract dofs definitions and dofs values
    dofs_ref = np.loadtxt(all_files[0], delimiter='=',
                          comments='      Variables:', usecols=0, dtype=str)
    nrs = len([r for r in dofs_ref if r[1] == 'R'])
    all_dofs = []
    for file in all_files:
        dofs = np.loadtxt(file, delimiter='=',
                          comments='      Variables:', usecols=0, dtype=str)
        assert (dofs == dofs_ref).all(), \
            f"{file} has different dofs than {all_files[0]}"
        dofs = np.loadtxt(file, delimiter='=',
                          comments='      Variables:', usecols=1)
        all_dofs.append(dofs)
    all_dofs = np.array(all_dofs) # [str][dof]

    # jump to furthest repeated structure and save the new order in new_set
    new_set = []
    j = 0
    while j < len(all_files) - 1:
        new_set.append(j)
        d_ij = (all_dofs[j] - all_dofs[j + 1:])
        # rescale distances and angles to have them in the same order of
        # magnitude. I could also evaluate the approximation for distances and
        # then for angles and use logicaland.
        condition = d_ij[:, nrs:] < -180
        while condition.any():
            d_ij[:, nrs:][condition] += 360
            condition = d_ij[:, nrs:] < -180
        condition = d_ij[:, nrs:] > 180
        while condition.any():
            d_ij[:, nrs:][condition] -= 360
            condition = d_ij[:, nrs:] > 180
        d_ij[:, nrs:] *= 1e-3 # trans 1 degree
        d_ij = abs(d_ij)
        close = np.isclose(d_ij, 0, atol=1e-3)
        by_struc = np.all(close, axis=1)
        same = np.where(by_struc)[0] + j
        if len(same) != 0:
            j = same[-1]
        j += 1
    new_set.append(len(all_files) - 1)

    # copy relevant files to a directory called subset
    all_files = np.array(all_files)[new_set]

    output_terminal("if [ ! -d subset ]; then mkdir subset; fi")
    for file in all_files:
        output_terminal("name=" + file + "; cp ${name%-dofs.dat}* subset")
    #output_terminal("cd subset; myutils rearange_force_files")
