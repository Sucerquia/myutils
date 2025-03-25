
==========
stretching
==========

.. container:: bash-script-title

   :ref:`[script] <stretching>` **myutils/sith/stretching.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool obtains the stretched configurations of a peptide by increasing the
      distance between carbons of the capping groups, constraining and optimizing
      using BMK exchange-correlation.
      
        -b  <number of breakages=1> The simulation will run until get this number of
            ruptures.
        -p  <peptide> One letter code of the amino acids forming the peptides. In
            this directory, a file called <peptide>-stretched00.pdb has to exist.
        -m  <method=0> index of stretching method. To see the options, use
            'myutils change_distance -h' to see the order.
        -r  restart stretching. In this case, this conde must be executed from
            the peptide's directory.
        -s  <size[A]=0.2> of the step that increases the distances.
      
        -v  verbose
        -h  prints this message.
      
