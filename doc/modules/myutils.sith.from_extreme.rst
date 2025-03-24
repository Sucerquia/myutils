myutils.sith.from\_extreme package
==================================

Submodules
----------


.. automodule:: myutils.sith.from_extreme.info_from_opt
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: myutils.sith.from_extreme
   :members:
   :undoc-members:
   :show-inheritance:

.. toctree::
   :hidden:

   bash_rsts_doc/opt_from_xyzs
   bash_rsts_scripts/opt_from_xyzs
   bash_rsts_doc/prepare_and_submit
   bash_rsts_scripts/prepare_and_submit
   bash_rsts_doc/extr_dofs
   bash_rsts_scripts/extr_dofs
   bash_rsts_doc/resubmit_failed
   bash_rsts_scripts/resubmit_failed
   bash_rsts_doc/opt_and_forces
   bash_rsts_scripts/opt_and_forces
   bash_rsts_doc/rearange_files
   bash_rsts_scripts/rearange_files
   bash_rsts_doc/workflow_from_extreme
   bash_rsts_scripts/workflow_from_extreme
   bash_rsts_doc/submit_forces_after_opt
   bash_rsts_scripts/submit_forces_after_opt
   bash_rsts_doc/reduce_ds
   bash_rsts_scripts/reduce_ds
   bash_rsts_doc/forces_from_xyzs
   bash_rsts_scripts/forces_from_xyzs
   bash_rsts_doc/after_optimization
   bash_rsts_scripts/after_optimization


after_optimization
------------------

.. mermaid::
   :align: center

   graph TD
   node1(["after_optimization"])
     node1 --> node2["Generate a path
                      of configurations
                      removing jumping pics
                      and interpolating if
                      the end_to_end distance
                      is larger than 0.2A.
                      **info_from_opt**"]
     click node2 "myutils.sith.from_extreme.html#myutils.sith.from_extreme.info_from_opt.info_from_opt" _self
     node2 --> node3["Extract dofs from all
                      xyz files using **extr_dofs**"]
     click node3 "myutils.sith.from_extreme.html#extr-dofs" _self
     node3 --> node4["Delete irrelevant changes and
                      adds intermedias to make DOFs
                      continuous with
                      **reduce_structs**"]
     click node4 "myutils.sith.from_extreme.html#myutils.sith.from_extreme.info_from_opt.reduce_structs" _self
     node4 --> node5["Creates a g09 com templete
                      to compute the forces of the
                      peptide based on \*-forces000.xyz
                      **forces_from_xyzs**"]
     click node5 "myutils.sith.from_extreme.html#forces-from-xyzs" _self
     node5 --> node6["Replaces the values of
                      DOFs in the template for
                      the continuous ones and
                      submits all calculations
                      with sbatch using
                      **single_g09**"]
     click node6 "myutils.bash_scripts.html#single-g09" _self


.. include:: bash_rsts_doc/after_optimization.rst

forces_from_xyzs
----------------

.. include:: bash_rsts_doc/forces_from_xyzs.rst

reduce_ds
---------

.. include:: bash_rsts_doc/reduce_ds.rst

submit_forces_after_opt
-----------------------

.. include:: bash_rsts_doc/submit_forces_after_opt.rst

workflow_from_extreme
---------------------

.. mermaid::
   :align: center

   graph TD
   node1(["workflow_from_extreme"])
     node1 --> node2["Find last xyz optimization
                      and start a new optimization
                      without constraints in a directory
                      called from_extreme"]
     node2 --> node3["used log output for
                      'after_optimization'"]
     click node3 "myutils.sith.from_extreme.html#after-optimization" _self

.. include:: bash_rsts_doc/workflow_from_extreme.rst

rearange_files
--------------

.. include:: bash_rsts_doc/rearange_files.rst

opt_and_forces
--------------

.. include:: bash_rsts_doc/opt_and_forces.rst

resubmit_failed
---------------

.. include:: bash_rsts_doc/resubmit_failed.rst

extr_dofs
---------

.. include:: bash_rsts_doc/extr_dofs.rst

prepare_and_submit
------------------

.. include:: bash_rsts_doc/prepare_and_submit.rst

opt_from_xyzs
-------------

.. include:: bash_rsts_doc/opt_from_xyzs.rst
