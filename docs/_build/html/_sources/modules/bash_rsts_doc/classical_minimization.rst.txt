
======================
classical_minimization
======================

.. container:: bash-script-title

   :ref:`[script] <classical_minimization>` **myutils/gromacs/classical_minimization.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool optimizes a configuration from an initial pdb file, which must be the
      first argument.
      
        -f  pdb file with the configuration to be minimized.
        -o  name of the output file with the optimized structure. Default same
            input (replaces the pdb of the input).
        -l  log file of the gromacs outputs. Default /dev/null
      
        -v  verbose.
        -h  prints this message.
      
