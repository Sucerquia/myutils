from myutils.g_values.g_valsetup import ext_xyz_from_npz, get_connectivity
import numpy as np
import glob
from ase import Atoms
from ase.io import write


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
        print(f'{npz_file}')
        data = np.load(npz_file)
        H = Atoms(symbols='H', positions=[data['hydrogen_xyz'][0]])
        radical = ext_xyz_from_npz(npz_file, create_xyz_files=False)[0]
        molecule = radical[:data['original_hydrogen_idxs'][0]] + H + \
            radical[data['original_hydrogen_idxs'][0]:]
        con = get_connectivity(molecule, 'vmd')
        if len(con[data['original_hydrogen_idxs'][0]]) != 1:
            write(f'totest.xyz', molecule)
            raise ValueError(f'Original hydrogen {data["original_hydrogen_idxs"][0]} ' +\
                f'should have only one bond. Got {len(con[data["original_hydrogen_idxs"][0]])}.')
            
        unique = []
        for i, connections in enumerate(con):
            for j in connections:
                if j > i:
                    unique.append([i, j])
        nb_ori = len(unique)

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
            assert len(unique) == nb_ori - 1, f'Number of bonds for original ' + \
                f"hydrogen {hi} is not correct. Expected {nb_ori - 1}, " +\
                f"got {len(unique)}. {len(unique), nb_ori - 1}"
            connectivity.append(unique)

        data = np.load(npz_file)
        data = {key: data[key] for key in data.files}
        data['bonds'] = np.array(connectivity)
        np.savez(npz_file, **data)

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
