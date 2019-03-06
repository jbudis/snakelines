Download reference
==================

At first, download sequences from NCBI according to selected genbank ids and prepare reference FASTA.
Also extract taxonomic ids and store them in the TAX file that is further used, particularly in metagenomics.


Purpose
-------

* Download genomic sequences and prepare reference database

Required inputs
---------------

* none

Generated outputs
-----------------

* Reference FASTA file with sequences

* Reference TAX file with taxonomies for FASTA file

Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/fasta_processing

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_download_reference.yaml

Example configuration:

.. literalinclude:: ../../example/fasta_processing/config_download_reference.yaml
   :language: yaml

Planned improvements
--------------------

* None


Included pipelines
------------------

* None