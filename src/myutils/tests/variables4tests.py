
from pathlib import Path
import numpy as np


gpa_endo_pdb = str(Path(__file__).parent) + "/references/GPA-endo.pdb"
gpa_exo_pdb = str(Path(__file__).parent) + "/references/GPA-exo.pdb"
gpa_endo_xyz = str(Path(__file__).parent) + "/references/GPA-endo.xyz"
gpa_opti_pdb = str(Path(__file__).parent) +  "/references/GPA-stretched00.pdb"
gpa_stre_pdb = str(Path(__file__).parent) +  "/references/GPA-stretched30.pdb"
frozendofs_dat = str(Path(__file__).parent) +  "/references/frozen_dofs.dat"

gpa_atoms = np.array(['CH3', '1HH3', '2HH3', '3HH3', 'C', 'O', 'N', 'H', 'CA',
                      'HA1', 'HA2', 'C', 'O', 'N', 'CD', 'HD1', 'HD2', 'CG',
                      'HG1', 'HG2', 'CB', 'HB1', 'HB2', 'CA', 'HA', 'C', 'O',
                      'N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C',
                      'O', 'N', 'H', 'CH3', '1HH3', '2HH3', '3HH3'])

gpa_residues = np.array(['ACE', 'ACE', 'ACE', 'ACE', 'ACE', 'ACE', 'GLY',
                         'GLY', 'GLY', 'GLY', 'GLY', 'GLY', 'GLY', 'PRO',
                         'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO',
                         'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'ALA',
                         'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA',
                         'ALA', 'ALA', 'NME', 'NME', 'NME', 'NME', 'NME',
                         'NME'])
