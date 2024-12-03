
.. container:: bash-script-title

   :ref:`[script] <resubmit_failed>` **myutils/sith/from_extreme/resubmit_failed.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Takes all the g09 log files given as arguments. If the log file does not report
      a proper termination, a new job is resubmitted (creating a backup first) using
      'myutils opt_and_forces'.
      
        -v  verbose
        -h  prints this message.
      
