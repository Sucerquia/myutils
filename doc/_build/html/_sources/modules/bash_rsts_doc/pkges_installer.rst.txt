
===============
pkges_installer
===============

.. container:: bash-script-title

   :ref:`[script] <pkges_installer>` **myutils/pkges_installer.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This tool installs all the next packages:
      
        -A  <path> ase*. Default=git repository when <path>='_'.
        -p  pymol.
        -v  vpython.
        -V  <path> vmol*. Default=git repository when <path>='_'.
        -n  ngl.
        -S  <path> sith*. Default=git repository when <path>='_'.
        -P  <path> pepgen*. Default=git repository when <path>='_'.
        -x  sphinx.
        -t  sphinx_rtd_theme.
      
        -h  prints this message.
      
      Consider the next options:
      
        -a  all packages.
        -d  <path>. directory where you want to store the packages installed
            from the source.
      
      Note: if you use the flag -a, all the packages are installed, even those with
      asterisk. The packages with asterisk are installed from the repository if you
      don't specify the path.
      
