
from pathlib import Path
import numpy as np


gpa_endo_pdb = str(Path(__file__).parent) + "/references/GPA-endo.pdb"
gpa_exo_pdb = str(Path(__file__).parent) + "/references/GPA-exo.pdb"
gpa_endo_xyz = str(Path(__file__).parent) + "/references/GPA-endo.xyz"
gpa_opti_pdb = str(Path(__file__).parent) +  "/references/GPA-stretched00.pdb"
gpa_stre_pdb = str(Path(__file__).parent) +  "/references/GPA-stretched30.pdb"
gpa_broken_xyz = str(Path(__file__).parent) +  "/references/GPA-broken.xyz"
frozendofs_dat = str(Path(__file__).parent) +  "/references/frozen_dofs.dat"

gpa_atoms = np.array(['CH3', '1HH3', '2HH3', '3HH3', 'C', 'O', 'N', 'H', 'CA',
                      'HA1', 'HA2', 'C', 'O', 'N', 'CD', 'HD1', 'HD2', 'CG',
                      'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'CA', 'HA', 'C', 'O',
                      'N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C',
                      'O', 'N', 'H', 'CH3', '1HH3', '2HH3', '3HH3'])

gpa_residues = np.array(['ACE ', 'ACE ', 'ACE ', 'ACE ', 'ACE ', 'ACE ', 'GLY ',
                         'GLY ', 'GLY ', 'GLY ', 'GLY ', 'GLY ', 'GLY ', 'PRO ',
                         'PRO ', 'PRO ', 'PRO ', 'PRO ', 'PRO ', 'PRO ', 'PRO ',
                         'PRO ', 'PRO ', 'PRO ', 'PRO ', 'PRO ', 'PRO ', 'ALA ',
                         'ALA ', 'ALA ', 'ALA ', 'ALA ', 'ALA ', 'ALA ', 'ALA ',
                         'ALA ', 'ALA ', 'NME ', 'NME ', 'NME ', 'NME ', 'NME ',
                         'NME '])

gpa_bonds = [[1, 2], [1, 3], [1, 4], [1, 5], [5, 6], [5, 7], [7, 8], [7, 9],
             [9, 10], [9, 11], [9, 12], [12, 13], [12, 14], [14, 15], [14, 24],
             [15, 16], [15, 17], [15, 18], [18, 19], [18, 20], [18, 21],
             [21, 22], [21, 23], [21, 24], [24, 25], [24, 26], [26, 27],
             [26, 28], [28, 29], [28, 30], [30, 31], [30, 32], [30, 36],
             [32, 33], [32, 34], [32, 35], [36, 37], [36, 38], [38, 39],
             [38, 40], [40, 41], [40, 42], [40, 43]]
