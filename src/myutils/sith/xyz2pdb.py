from ase.io import read
import glob


# add2executable
def xyz2pdb(xyzfile, pdbtemplate, pdbfile=None, withx=0):
    if pdbfile is None:
        pdbfile = xyzfile.split('.')[0] + '.pdb'
    atoms = read(xyzfile)
    positions = atoms.positions
    with open(pdbtemplate) as f:
        lines = f.readlines()

    i_atom = 0
    for index in range(len(lines)):
        line_split = lines[index].split()
        if 'ATOM' in line_split[0]:
            atomposition = positions[i_atom]
            i_atom += 1
            for i in [0, 1, 2]:
                toreplace = line_split[6+withx+i]
                n2dot_toreplace = len(toreplace.split('.')[0])
                replacefor = "{:.3f}".format(atomposition[i])
                n2dot_replacefor = len(replacefor.split('.')[0])

                if n2dot_replacefor > n2dot_toreplace:
                    toreplace = ' '*(n2dot_replacefor - n2dot_toreplace) + \
                        toreplace
                if n2dot_toreplace > n2dot_replacefor:
                    replacefor = ' '*(n2dot_toreplace - n2dot_replacefor) + \
                        replacefor
                lines[index] = lines[index].replace(toreplace, replacefor)

    with open(pdbfile, "w") as f:
        [f.write(line) for line in lines]


# add2executable
def all_xyz2pdb(template):
    configs = glob.glob('*.xyz')
    configs.sort()
    configs.pop(0)
    for config in configs:
        transform(config, template)

    return 0
