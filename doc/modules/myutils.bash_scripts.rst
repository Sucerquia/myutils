
.. container:: bash-script-title

   **myutils/bash_scripts/single_g09.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This code runs one optimization using gaussian in one of the clusters. You
      have to create the input file and give it (without .com extension) as first
      argument when run this code.
      
        -f  <com file> name if the gaussian input file without extension (.com).
        -c  run in server.
      
        -h  prints this message.
      

.. container:: bash-script-title

   **myutils/bash_scripts/bash-template.sh**

.. container:: bash-script-doc

   .. line-block::
      
      Use this template to create your scripts with a standard structure
      
        -d  <variable> add the description of the variable.
        -c  usually used when submitted in a cluster to import modules.
      
        -v  verbose.
        -h  prints this message.
      

.. container:: bash-script-title

   **myutils/bash_scripts/find_blocks.sh**

.. container:: bash-script-doc

   .. line-block::
      
      This code extracts the sections in a file starting and finishing with specific
      patterns without including the lines containing those patterns. Check the next
      options:
      
        -f  <file> file that shows 
        -s  <pattern> pattern that defines the beginning of the block. This line
            is not included in the block.
        -e  <pattern> pattern that defines the end of the block. This line is not
            included in the block.
        -i  use this flag if the start and the end are indexes
        -o  <output='output'> 'terminal' or the name of the output without
            extension. In the later case, the output will be stored in a file
            called <output>.dat if the flag -i is given or in files called
            <output>_<n>.dat, where n is the number of appearence of the block in
            the file.
      
        -v  verbose of what's the code doing.
        -h  prints this message.
      
