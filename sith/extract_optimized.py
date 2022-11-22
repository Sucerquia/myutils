from ase.io import read, write
from glob import glob

files = glob('*.log')
files.sort()
trajectory = []

for opt in files:
    config = read(opt)
    write('opt-'+opt.split('.')[0]+'.pdb', config)
    trajectory += [config]

write('opt-all.pdb', trajectory)   
    
