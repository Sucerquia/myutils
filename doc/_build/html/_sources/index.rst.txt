.. _myutils:

=======
MyUtils
=======

MyUtils is a set of tools that you can use for scientific proporsals, mainly
related with the field of simulations of molecules, although there are some
parts that could be interesting for all the scientific community. This set of
tools was created and it is currently maintained by Daniel Sucerquia, PhD
student at the Heidelberg Institute for Theoretical Studies (HITS) and the
Max Planck Institute for Polymer Research (MPIP).

This is a general overview of the kind of tools you can find in this project:

.. mermaid::
   :align: center

   graph TD
   node1["myutils"]
   click node1 "modules/myutils.html" _self
     node1 --> node2["ase_utils"]
     click node2 "modules/myutils.ase_utils.html" _self
     node1 --> node3["bash_scripts"]
     click node3 "modules/myutils.bash_scripts.html" _self
     node1 --> node4["cli"]
     click node4 "modules/myutils.cli.html" _self
       node4 --> node5["pkg_structure"]
       click node5 "modules/myutils.cli.pkg_structure.html" _self
     node1 --> node6["gromacs"]
     click node6 "modules/myutils.gromacs.html" _self
     node1 --> node7["pre-deprected"]
     click node7 "modules/myutils.pre-deprected.html" _self
     node1 --> node8["sith"]
     click node8 "modules/myutils.sith.html" _self
       node8 --> node9["from_extreme"]
       click node9 "modules/myutils.sith.from_extreme.html" _self

-------
Content
-------

.. toctree::
   about
   install
   modules/myutils
   tutorials/tutorials



