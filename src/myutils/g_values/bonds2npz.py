from myutils.g_values.g_valsetup import ext_xyz_from_npz, get_connectivity
import numpy as np
import glob
from ase import Atoms
from ase.io import write

# add2executable
def addfile2npz(npz_file_ref, npz_file_add, file_name):
    """
    Add file information to npz files from another one

    Parameters
    ==========
    npz_file_ref: str
        Path to the reference npz file.
    npz_file_add: str
        Path to the npz file to add information.

    Returns
    =======
    (dict) data in npz file.
    """
    ref = np.load(npz_file_ref)
    add = np.load(npz_file_add)
    add_data = {key: add[key] for key in add.files}
    add_data[file_name] = ref[file_name]
    np.savez(npz_file_add, **add_data)
    return add_data


# add2executable
def bonds2npz(npz_file):
    """
    Add bond information to npz files. The bond information is obtained from
    the xyz files that are created from the npz files. The connectivity is
    obtained using the VMD method.
    
    Parameters
    ==========
    npz_file : str
        Path to the npz file or directory containing npz files. It could also
        be a path to a .dat file that contains a list of npz files.

    Returns
    =======
    (list) connectivity information for each molecule in the npz file(s).
    """
    if npz_file.endswith('.npz'):
        npz_files = [npz_file]
    elif npz_file.endswith('.dat'):
        npz_files = np.loadtxt(npz_file, dtype=str)
    else:
        npz_files = glob.glob(npz_file + '/*.npz')

    for npz_file in npz_files:
        try:
            print(f'{npz_file}')
            data = np.load(npz_file)
            molecule = Atoms(numbers=data['original_atomic_numbers'],
                            positions=data['original_xyz'])

            con = get_connectivity(molecule, 'vmd')

            # hydrogen should have only one bond.
            if len(con[data['original_hydrogen_idxs'][0]]) != 1:
                raise ValueError(f'Original hydrogen {data["original_hydrogen_idxs"][0]} ' +\
                    f'should have only one bond. Got {len(con[data["original_hydrogen_idxs"][0]])}.')

            connectivity = []
            for hi in data['original_hydrogen_idxs']:
                unique = []
                for i, connections in enumerate(con):
                    if i == hi:
                        continue
                    elif i < hi:
                        l = i
                    else:
                        l = i - 1
                    for j in connections:
                        if j == hi:
                            continue
                        elif j < hi:
                            m = j
                        else:
                            m = j - 1
                        if j > i:
                            unique.append([l, m])
                connectivity.append(unique)
            data = np.load(npz_file)
            data = {key: data[key] for key in data.files}
            data['bonds'] = np.array(connectivity)
            np.savez(npz_file, **data)
        except:
            with open('failed.dat', 'a') as f:
                f.write(f'{npz_file}\n')

    return connectivity


# add2executable
def molid2npz(npz_file):
    """
    Add molecular ID information to npz files. The molecular ID information
    
    Parameters
    ==========
    npz_file : str
        Path to the npz file or directory containing npz files.
    
    Returns
    =======
    
    """
    if npz_file.endswith('.npz'):
        npz_files = [npz_file]
    else:
        npz_files = glob.glob(npz_file + '/*.npz')
        npz_file = [f[2:] for f in npz_files]
    
    for npz_file in npz_files:
        data = np.load(npz_file)
        data = {key: data[key] for key in data.files}
        data['molid'] = np.array([f'{npz_file[:-4]}_{i:03d}' for i in range(len(data['xyz']))])
        np.savez(npz_file, **data)
