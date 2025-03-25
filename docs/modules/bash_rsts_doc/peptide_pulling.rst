
===============
peptide_pulling
===============

.. container:: bash-script-title

   :ref:`[script] <peptide_pulling>` **myutils/gromacs/peptide_pulling.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool creates the trajectory of a given peptide pulled by an external force.
      Consider the next options:
      
        -a  properties you want to analyse. For example "-d -r". Default "-d -L".
            For more information, check: myutils analysis -h
        -f  forces to stretch the peptide in [kJ mol^-1 nm^-1]. eg 100,200.
            Default 200
        -g  gromacs binary. For example gmx or gmx_mpi. Default gmx.
        -o  pepgen flags
        -p  peptide.
      
        -v  verbose.
        -h  prints this message.
      
