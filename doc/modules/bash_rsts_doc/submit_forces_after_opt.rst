
.. container:: bash-script-title

   :ref:`[script] <submit_forces_after_opt>` **myutils/sith/from_extreme/submit_forces_after_opt.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Searches all the <dirs>/forces/*-opt.log files where <dirs> are all the
      directories in the current location. With those files a job called
      <file_name>_forces is submitted with sbatch using 'myutils compute_forces'
      
        -v  verbose
        -h  prints this message.
      
