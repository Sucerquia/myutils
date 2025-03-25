
==============
doc_pythonfile
==============

.. container:: bash-script-title

   :ref:`[script] <doc_pythonfile>` **myutils/cli/pkg_structure/doc_pythonfile.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Code that automatically creates or update the '.rst' file for the documentation
      of all python classes and functions that finds in a python file. Also add the
      file to <path_to_doc>/modules.rst if it does not exist.
      
        -f  <file> file name with the stem relative to the directory that contains
            all the python codes to be documented (pkg_path, see -p below).
        -m  <mod_doc_path> path to the directory that stores the modules
            documentation.
        -n  <pkg_name> name of the package to be documented.
        -p  <pkg_path> path to directory that contains all the python codes to be
            documented. Usually src directory.
        -h  prints this message.
      
      Note: This documentation is used in doc_modules.sh
      
