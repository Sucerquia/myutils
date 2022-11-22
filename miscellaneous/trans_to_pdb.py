''' 
This code takes the last configuration of one file (e. one gaussian *.log file)
and saves the last configuration in a pdb file.

python trans_G09-pdb.py arg1 <arg2>

arg1: file that contains the last configuration.

arg2: (optional) name of the output without expension.

Note: In case of unexisting argument 2, the output will have the same name
that argument 1 but with extension .pdb

example:
python trans_G09-pdb.py optimization.log
'''

from ase.io import read, write
import glob
import sys


def change_format(names):
    """"
    This function takes the last configuration of one file (e. one gaussian 
    *.log file) and saves the last configuration in a pdb file.

    Parameters
    ==========
    names is a one/two-component list, the first element must be the string of
    the files to be modified (e.g. './*.log' ). The second element (optional) 
    contains the output name. In case to receive one-component element, this
    function will save the new file with the same name but pdb extension
    """

    to_modify = glob.glob(names[0])
    to_modify.sort()
    n_files_to_modify = len(to_modify)


    if len(names) > 1:  # there is an output name
        if n_files_to_modify > 1:  # there are several files to be modified
            output = [names[1]+str(i)+'.pdb' for i in range(len(to_modify))]
        else:
            output = [names[1]+'.pdb']
    else:
        output = list()
        for name in to_modify:
            index_rename = name.rfind('.')
            rename = name[:index_rename]
            output.append(rename+'.pdb')

    for i in range(n_files_to_modify):
        print(f" {to_modify[i]} ---> {output[i]} ")
        a = read(to_modify[i], index=-1)
        write(output[i], a)
        

if __name__ == '__main__':
    change_format(sys.argv[1:])