
.. container:: bash-script-title

   :ref:`[script] <add_python_doc>` **myutils/cli/pkg_structure/add_python_doc.sh**

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
      
