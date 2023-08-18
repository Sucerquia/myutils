# Intro
set of tools for different aims. Created by Daniel Sucerquia, PhD student at HITS, 2023.

Twitter (X): @sucer6

# Requirements
- python >= 3.9
- pytest
- ASE
- pymol
- pepgen
- ngl_view

Optional:
- sphinx
- sphinx_rtd_theme


# Installation with conda

```bash
conda create --name myutils python=3.9
conda activate myutils
pip install git+https://gitlab.com/ase/ase.git
conda install -c conda-forge pymol-open-source
pip install git+https://github.com/hits-mbm-dev/pepgen.git

pip install git+https://github.com/Sucerquia/myutils.git
```

Note: all the codes that are made to run in a cluster asume that the
environment "myutils" exist and has this package already installed.


# Tests
2COMPLETE define a direct way to run the internal tests.

There are some markers for specific cases, to activate them, run the test as:

2COMPLETE after defining a way to run the tests directly from the package, rewrite this
pytest -m "marker1 or marker2 or not marker3 ..."

where the markers could be one of the next:

-  runincluster :  some tests are done to submit jobs in the queue sysrem
slrum. You can run these tests in an interactive session using 2COMPLETE
-  toolong : tests that last more than one minute running. They are usually
codes with g09 processes running in 8 cores.