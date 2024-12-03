
.. container:: bash-script-title

   :ref:`[script] <workflow_from_extreme>` **myutils/sith/from_extreme/workflow_from_extreme.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool creates the files to do the sith analysis by optimizing a molecule
      that was just about to get a first rupture, then takes the intermedia steps and
      find the internal forces. Consider the next options:
      
        -c  run in cascade. (modules are loaded)
        -p  <peptide>. directory or xyzfile of last conf. Chains of aminoacids to
            be evaluated. For example, "./AAA/" would optimize the last
            stretched a trialanine peptide (where last means after organizing
            alphabetically).
        -l  <number of amino acids in the peptide> It will be assumed that the
            xyz file starts with the letter code of the amino acids.
      
        -v  verbose.
        -h  prints this message.
      
      Note: it is assumed that the file of the last configuration is named as:
      <amino acids-code>-<description><number of stretching>.xyz
      
