
========
analysis
========

.. container:: bash-script-title

   :ref:`[script] <analysis>` **myutils/gromacs/analysis.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool helps you to extract information from a gromacs trajectory. You have
      to give the name of the file, the gromacs binary and the flag of the function
      you want to compute considering the next options:
      
        -a  computes all the properties available in this tool and merges them.
        -c  extract the last configuration of the trajectory.
        -d  computes distance between the atoms specified 'distance' in the index
            file.
        -e  computes potential energy.
        -f  name of the file without extension. edr, trr and gro file have to have
            the same name.
        -g  gromacs binary. usualy gmx.
        -k  keeps the copies of the files with the same name. Default: false.
        -L  extracts the largest configuration.
        -l  log file of the gromacs outputs.
        -m  merges all data in one file.
        -r  computes ramachandran angles.
        -s  computes the energy of the subsystem in the trajectory, in this case,
            the protein.
      
        -v  verbose.
        -h  prints this message.
      
