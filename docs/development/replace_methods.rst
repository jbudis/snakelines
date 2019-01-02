Replace default methods
=======================

Methods in each analysed step can be easily replaced with alternative tool or custom solution.
All supported values are enumerated in the configuration value.
For example, you may use spades or unicycler for assembly by changing value under `method:`

.. code-block:: yaml

   assembly:
       assembler:                      # Method for joining reads
           method: spades              # Supported values: spades, unicycler
           mode: standard              # Supported values: standard, meta, plasmid, rna, iontorrent
           careful: True               # Tries to reduce number of mismatches and short indels, longer runtime


You may use tool that is not yet implemented, or custom solution.
In that case, you need to create new rule in the /rules/assembly/assembler/<your_assembler>.snake.
Rule should contain the same set of outputs as the other assemblers as is defined in the dependency file (/src/dependency.yaml).
Than you can easily switch analysis to your solution changing configuration file:

.. code-block:: yaml

   assembly:
       assembler:                       # Method for joining reads
           method: <your_assembler>     # Supported values: spades, unicycler, <your_assembler>
           <param_name>: <its_value>    # Specific attributes for <your_assembler>

Do not forget to change configuration of the tool according to your parameters as well.
