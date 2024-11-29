
.. container:: bash-script-title

   **myutils/cli/pkg_structure/add_python_doc.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Code that explores the files in the package and automatically creates the
      documentation of all classes and functions that finds in it.
      
        -d  <dir1,dir2...> directories to be ignored.
            Default: 'pycache,tests,tutorials,pre-deprected'
        -f  <fil1,fil2...> files to be ignored. Default: '__init__'
        -n  <name> pkg name. Default: myutils
        -p  <absolute_path> path directory to be checked (no relative path).
            Default: "$myutils path"
      
        -v  verbose.
        -h  prints this message.
      

.. container:: bash-script-title

   **myutils/cli/pkg_structure/bash_style.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Check bash style of all bash files in a directory, and its subdirectories.
      
        -d  directory. Default: "$myutils -path"
      
        -v  verbose.
        -h  prints this message.
      

.. container:: bash-script-title

   **myutils/cli/pkg_structure/doc_modules.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Code that explores the files in the package and automatically create the
      documentation of all classes and functions that finds in it.
      
         -d   <dir1,dir2...> directories to be ignored. Default: 'tests,cli'
         -f   <fil1,fil2...> files to be ignored. Default: '__init__'
         -p   <absolute_path> path directory to be checked (no relative path).
              Default: "$myutils path"
         -m   <mod_doc_path> absolute path to the directory that stores the modules
              documentation. Default: <mod_path>/../../doc/modules
         -n   <name> pakage name. Default: myutils
      
         -h   prints this message.
      

.. container:: bash-script-title

   **myutils/cli/pkg_structure/check_structure.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Check the structure of a package. All checkers run by default.
      
        -d  src directory of the package. Defatul: "$myutils -path"
        -p  pep8 convention in all python scripts.
        -s  ShellCheck in all bash scripts.
        -t  check tests.
      
        -h  prints this message.
      

.. container:: bash-script-title

   **myutils/cli/pkg_structure/python_doc_fixer.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Code takes the documentation of a function and checks the documentation adding
      TODOs in the missing parts. The output is stored in a file called
      final_<method>_doc.txt
      
        -f  <method> Function to be checked
        -m  <module> Module that contains the Function
        -s  <n_spaces> number of leading spaces.
      
        -v  verbose.
        -h  prints this message.
      

.. container:: bash-script-title

   **myutils/cli/pkg_structure/doc_pythonfile.sh**

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
      

.. container:: bash-script-title

   **myutils/cli/pkg_structure/bash_basic_structure.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Check the basic structure of all the given files as arguments.
      
      It checks that all the bash scripts have flags support, verbose and help.
      
        -v  verbose.
        -c  run in a cluster.
      
