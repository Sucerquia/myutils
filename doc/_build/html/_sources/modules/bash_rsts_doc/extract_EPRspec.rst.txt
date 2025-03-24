
===============
extract_EPRspec
===============

.. container:: bash-script-title

   :ref:`[script] <extract_EPRspec>` **myutils/g_values/analysis/extract_EPRspec.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Extract the EPR absorption spectrum from an orca output file.
      
        -e  <file.dat> experimental field vs absorption file.
        -f  Use this flag to take into account hyperfine corrections. This uses a lot
            of RAM memory. Be sure that you have enough memory or that you filtered
            the nuclei to compute hyperfine correction.
        -O  <file.out='EPRII_i.o'> orca output file with computed EPR quantities. 
        -o  <output.dat='spectrum_wo_hyFiCorr.dat'> dat output file where you want to
            save the field vs spectrum.
        -m  <float=179.813> experimental value of the microwave frequency. The
            default value corresponds to G-band experiments.
        -n  <int=501> number of data points used to predict the absorption spectrum. 
      
        -v  verbose.
        -h  prints this message.
      
