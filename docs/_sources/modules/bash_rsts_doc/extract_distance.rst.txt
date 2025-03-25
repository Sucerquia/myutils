
================
extract_distance
================

.. container:: bash-script-title

   :ref:`[script] <extract_distance>` **myutils/gromacs/extract_distance.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Extract distance between two atoms from gromacs trajectory.
      
        -r  <res1,res2>, residues indexes of the atoms to compute the distance.
        -a  <a1, a2>, names of the atoms to compute the distance.
        -t  <traj file> path to the trajectory file to extract the distance.
        -o  <output> name of the output without extension (.dat)
        -g  <gro_file> gro file used to create the trajectory.
      
        -v  verbose.
        -h  prints this message.
      
