
============
gval_postpro
============

.. container:: bash-script-title

   :ref:`[script] <gval_postpro>` **myutils/g_values/analysis/gval_postpro.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Create all the files after complete the computation of the gvalues with Orca.
      
        -c  <molecule candidate> Name of the radical candidate. It is assumed that
            the files related with this candidate are in a directory with the same
            name.
        -E  <experimental values> guess of the gvalues obtained experimentally in
            python list format, f.e. '[2.0062, 2.0055, 2.0022]'
        -e  <file.dat> experimental field vs absorption file.
        -f  Use this flag to take into account hyperfine corrections. This uses a lot
            of RAM memory. Be sure that you have enough memory or that you filtered.
        -O  <file.out='EPRII_i.o'> orca output file with computed EPR quantities.
        -o  <output.dat='spectrum_wo_hyFiCorr.dat'> dat output file where you want to
            save the field vs spectrum.
        -m  <float=179.813> experimental value of the microwave frequency. The
            default value corresponds to G-band experiments.
        -n  <int=501> number of data points used to predict the absorption spectrum
      
        -v  verbose.
        -h  prints this message.
      
      This code should produce:
      
       - A render of each model called opt_EPRII.png
       - The computed spectrum obtained by easyspin in a file called spectrum_wo_hyFiCorr.dat
       - A plot of the gvalues of all the candidates called gvalues.png
       - A table of the gvalues in gvalues_table.md
       - A file called vmd_image.md in each one of the directories of the candidates.
      
      This code should be executed in the folder containing the candidate.
      
