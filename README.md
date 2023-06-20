# Intro
set of tools for different aims. Created by Daniel Sucerquia, PhD student at HITS.

Twittwe: @sucer6

# Notes for developers
definition:2COMPLETE is a marker for the developers to know what is incomplete. please add
those parts also as en issue in the github repository.

# Requirements
- python >= 3.9
- Numpy
- pytest
- ASE
- pymol and pepgen (optional: only for sith sources codes)
- ngl_view (optional: only for molecules visualization.)
- mathplotlib (optional: only for plottlers and molecules visualization)

# Installation with conda
conda create --name myutils python=3.9
conda activate myutils
pip install -e .
conda install -c conda-forge pymol-open-source
pip install -e .

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