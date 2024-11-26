import numpy as np


def index_dof(sith, target: tuple):
    """
    Search the index of a specific dof.

    Parameter
    =========
    target: np.ndarray
        degree of freedom.

    Return
    ======
    (int) index
    """
    for i, dof in enumerate(sith.dim_indices):
        dof_wo_0 = dof[np.nonzero(dof)[0]]
        target_wo_0 = target[np.nonzero(target)[0]]

        if len(dof_wo_0) != len(target_wo_0):
            continue

        if (dof_wo_0 == target_wo_0).all() or \
            (dof_wo_0 == target_wo_0[::-1]).all():
            return i
    raise ValueError("Non-found dof.")


def organize_lengths(sith1, sith2):
    n_distances = sith2.dims[1]

    reorder = []
    for dof in np.array(sith1.dim_indices[:n_distances]):
        new_index = index_dof(sith2, dof)
        reorder.append(new_index)

    return np.array(reorder)


def extract_common(sith1, sith2):
    """
    Finds the DOfs that are in common in both sith objects

    Return
    ======
    (tuple), list of indexes in sith1 that are in sith2, organized list 
    """
    sith1_indx = []
    sith2_indx = []
    for i, dof in enumerate(np.array(sith1.dim_indices)):
        try:
            new_index = index_dof(sith2, dof)
            sith1_indx.append(i)
            sith2_indx.append(new_index)
        except ValueError:
            continue

    return np.array(sith1_indx), np.array(sith2_indx)
