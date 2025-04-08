.. _install:

=======
Install
=======

------------
Dependencies
------------

For different parts of the package, different tools are required. None of the
packages listed below are really mandatory to use Myutils. You can just install
MyUtils and install the required dependencies once you get an error like modules
not found. This is just a heads up that this is a very likely error for
first-time users of MyUtils.

------------
Installation
------------

Run these lines on your terminal:

.. code-block:: bash

    git clone git@github.com:Sucerquia/myutils.git
    cd myutils
    pip install -e .

Note that "-e" flag allows to use MyUtils in such a way that you can make
changes on it and use them immediately without having to reinstall the package
again. However, you can avoid this flag and a copy of the package would be
generated in your python packages path.

.. _orca-setup:

-----------
Orca set up
-----------

To use the tools of MyUtils that include the use of orca in any of the steps, you
have to create a directory file called :code:`$mine/sw/orca/setup_orca.sh` containing
at least the next lines

.. code-block:: bash

    orca_path="<absolute path to the Orca directory>"
    export PATH=${orca_path}:$PATH
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$orca_path"
    
    orca="$orca_path/orca"

Note that it also applies if you run orca in a cluster. In that case or if it is
done by default (in your bash script for example), you only have to define "orca"
variable as the binary of orca (last line in the script above).