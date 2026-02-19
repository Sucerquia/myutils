from myutils.g_values.g_valsetup import ext_xyz_from_npz, get_connectivity
import numpy as np
import glob

# add2executable
def bonds2npz(npz_file):
    """
    Add bond information to npz files. The bond information is obtained from
    the xyz files that are created from the npz files. The connectivity is
    obtained using the VMD method.
    
    Parameters
    ==========
    npz_file : str
        Path to the npz file or directory containing npz files.
    
    Returns
    =======
    (list) connectivity information for each molecule in the npz file(s).
    """
    if npz_file.endswith('.npz'):
        npz_files = [npz_file]
    else:
        npz_files = glob.glob(npz_file + '/*.npz')
        npz_files = [f[2:] for f in npz_files]

    for npz_file in npz_files:
        molecules = ext_xyz_from_npz(npz_file, create_xyz_files=False)
        connectivity = []
        for atoms in molecules:
            con = get_connectivity(atoms, 'vmd')
            unique = []
            for i, connections in enumerate(con):
                for j in connections:
                    if j > i:
                        unique.append([i, j])
            connectivity.append(unique)
        

        data = np.load(npz_file)
        data = {key: data[key] for key in data.files}
        data['bonds'] = np.array(connectivity)
        np.savez(npz_file, **data)

    return connectivity
# %%

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
        npz_files = [f[2:] for f in npz_files]
    
    for npz_file in npz_files:
        data = np.load(npz_file)
        data = {key: data[key] for key in data.files}
        data['molid'] = np.array([f'{npz_file[:-4]}_{i:03d}' for i in range(len(data['xyz']))])
        np.savez(npz_file, **data)
