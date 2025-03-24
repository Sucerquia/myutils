myutils.sith package
====================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   myutils.sith.from_extreme

Submodules
----------


.. automodule:: myutils.sith.analysis
   :members:
   :undoc-members:
   :show-inheritance:


.. automodule:: myutils.sith.compare_siths
   :members:
   :undoc-members:
   :show-inheritance:


.. automodule:: myutils.sith.g09_xyz
   :members:
   :undoc-members:
   :show-inheritance:


.. automodule:: myutils.sith.protonate
   :members:
   :undoc-members:
   :show-inheritance:


.. automodule:: myutils.sith.sith_plots
   :members:
   :undoc-members:
   :show-inheritance:


.. automodule:: myutils.sith.sith_tools
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: myutils.sith
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :hidden:

   bash_rsts_doc/swap_atoms_in_com
   bash_rsts_scripts/swap_atoms_in_com
   bash_rsts_doc/workflow
   bash_rsts_scripts/workflow
   bash_rsts_doc/compute_forces
   bash_rsts_scripts/compute_forces
   bash_rsts_doc/find_forces
   bash_rsts_scripts/find_forces
   bash_rsts_doc/clean_ds
   bash_rsts_scripts/clean_ds
   bash_rsts_doc/proline_mod
   bash_rsts_scripts/proline_mod
   bash_rsts_doc/extract_forces
   bash_rsts_scripts/extract_forces
   bash_rsts_doc/stretching
   bash_rsts_scripts/stretching


stretching
----------

.. mermaid::
   :align: center

   graph TD
   node1(["stretching"])
     node1 --> node2["Optimization of the input
                      structure"]
     node2 --> node3["Takes last optimization
                      stretches it with
                      'change_distance'"]
     click node3 "myutils.ase_utils.html#myutils.ase_utils.tools.change_distance" _self
     node3 --> node4["Add preferences to g09
                      input and run it. it's
                      checked if the
                      convergence worked; it's
                      rerun in case of one non
                      convergence"]
     node4 --> node5["Extract xyz from the
                      optimization process with
                      'log2xyz'"]
     click node5 "myutils.sith.html#myutils.sith.g09_xyz.log2xyz" _self
     node5 --> node6["Extract diffrent bonds in
                      respect to the previous
                      optimization with
                      'diff_bonds'"]
     click node6 "myutils.ase_utils.html#myutils.ase_utils.tools.diff_bonds" _self
     node6 --> node7{"did
                      a bond
                      dissapear?"}
     node7 --> |No| node3
     node7 --> |Yes| node8{"enough
                            ruptures?"}
     node8 --> |No| node9["save broken structure to
                           rupture directory with
                           'create_bck'."]
     click node6 "myutils.ase_utils.html#myutils.ase_utils.tools.diff_bonds" _self
     node9 --> node3
     node8 --> |Yes| node10(["finish"])

.. include:: bash_rsts_doc/stretching.rst

extract_forces
--------------

.. include:: bash_rsts_doc/extract_forces.rst

proline_mod
-----------

.. include:: bash_rsts_doc/proline_mod.rst

clean_ds
--------

.. include:: bash_rsts_doc/clean_ds.rst

find_forces
-----------

.. include:: bash_rsts_doc/find_forces.rst

compute_forces
--------------

.. include:: bash_rsts_doc/compute_forces.rst

workflow
--------

.. mermaid::
   :align: center

   graph TD
   node1(["Workflow"])
     node1 ==> node2["Creates the peptide with pepgen and minimizes it with
                      'classical_minimization'"]
     click node2 "myutils.gromacs.html#classical-minimization" _self
     node2 --> node3["'proline_mod'"]
     click node3 "myutils.sith.html#proline-mod" _self
     node3 --> node4["'protonate'"]
     click node4 "myutils.sith.html#myutils.sith.protonate.protonate" _self
     node4 --> node5["**'stretching'**"]
     click node5 "myutils.sith.html#stretching" _self
     node1 -.->|restart| node5
     node5 --> node6["'classical_energies'"]
     click node6 "myutils.gromacs.html#classical-energies" _self
     node5 --> node7["'find_forces'"]
     click node7 "myutils.sith.html#find-forces" _self
     node5 --> node8["**'workflow_from_extreme'**"]
     click node8 "myutils.sith.from_extreme.html#workflow-from-extreme" _self

.. include:: bash_rsts_doc/workflow.rst

swap_atoms_in_com
-----------------

.. include:: bash_rsts_doc/swap_atoms_in_com.rst
