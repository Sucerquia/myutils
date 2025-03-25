
==============
extract_forces
==============

.. container:: bash-script-title

   :ref:`[script] <extract_forces>` **myutils/sith/extract_forces.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Extract the forces and indexes of the DOFs from the log files (g09). The output
      is a set of files called <pep>-forces<n_stretching>.fchk containing the
      information in fchk g09 format.
      
        -d  <path>. directory where forces_files.log are located. Default ./forces
      
        -v  verbose.
        -h   prints this message.
      
