
=======
pulling
=======

.. container:: bash-script-title

   :ref:`[script] <pulling>` **myutils/gromacs/pulling.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This code executes the pulling adding an external force along the x axis to the
      carbon atoms of the NME and ACE caps. Consider the next options:
      
        -f  forces to stretch the peptide in [kJ mol^-1 nm^-1].
        -g  gromacs binary. For example gmx or gmx_mpi. Default gmx.
        -l  log file of the gromacs outputs. Set it as /dev/tty to print the gromacs
            output on the terminal. Default /dev/null
        -s  steps in the MD pulling. Default: 10000
      
        -v  verbose.
        -h  prints this message.
      
      Note that the file equilibration must exist.
      
