Infer phylogeny between genomic sequences
=========================================

Compare sequences from the reference FASTA file and visualise their relationship in form of phylogenetic tree and
interactive multiple alignment.


Purpose
-------

* Assess relationship and similarity between genomic sequences

Required inputs
---------------

* Reference genome in fasta format

::

   |-- reference/<reference>
           |-- <reference>.fa

Generated outputs
-----------------

* Phylogenetic tree with distances between genomic sequences

* Interactive visualization of multiple alignment of reference sequences

Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/fasta_processing

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_infer_phylogeny.yaml \
      -- use-conda

Example configuration:

.. literalinclude:: ../../example/fasta_processing/config_infer_phylogeny.yaml
   :language: yaml

Planned improvements
--------------------

* SVG figure size should be scaled with the number of the sequences


Included pipelines
------------------

* None