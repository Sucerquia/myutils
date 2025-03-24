
===========
doc_modules
===========

.. container:: bash-script-title

   :ref:`[script] <doc_modules>` **myutils/cli/pkg_structure/doc_modules.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Code that explores the files in the package and automatically create the
      documentation of all classes and functions that finds in it.
      
         -d   <dir1,dir2...> directories to be ignored. Default: 'tests,cli'
         -f   <fil1,fil2...> files to be ignored. Default: '__init__'
         -p   <relative_path> relative path directory to be checked. Relative in
              respect to the directory that stores the modules documentation (see the
              flag -m). Default=../../src/myutils
         -m   <absolute_path> absolute path to the directory that stores the modules
              documentation. Default: $(myutils path)/../../doc/modules
         -n   <name> pakage name. Default: myutils
      
         -h   prints this message.
      
