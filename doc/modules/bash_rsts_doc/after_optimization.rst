
.. container:: bash-script-title

   :ref:`[script] <after_optimization>` **myutils/sith/from_extreme/after_optimization.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Creates the com files from the xyz structures extracted from a g09 log file and
      submit the corresponding jobs to compute the forces.
      
        -l  <log_file> optimization g09 logfile.
        -n  <name> standard name. Usually pep name.
      
        -h  prints this message.
      
