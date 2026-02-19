# %%
from myutils.g_values.g_valsetup import ext_xyz_from_npz, get_connectivity
import numpy as np

# add2executable
def bonds2npz(npz_file):
    data = np.load(npz_file)
    molecules = ext_xyz_from_npz(npz_file, create_xyz_files=False)
    connectivity = []
    for atoms in molecules:
        con = get_connectivity(atoms)
        connectivity.append(con)

    return connectivity
# %%