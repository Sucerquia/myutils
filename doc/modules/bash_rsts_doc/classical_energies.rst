
.. container:: bash-script-title

   :ref:`[script] <classical_energies>` **myutils/gromacs/classical_energies.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Tool that computes the classical energy from a set of pdb files in the
      executing directory.
      
        -l  log file of the gromacs outputs. Default /dev/null
        -n  use this flag to NOT transform all xyz files into pdbs. In this case is
              assumed that the pdbs already exist.
      
        -v  verbose.
        -h  prints this message.
      
