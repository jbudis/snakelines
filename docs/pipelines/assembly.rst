Assemble reads
==============

Join reads with overlaps into larger continuous genomic sequences, contigs.


Purpose
-------

* Determine genomic sequence of sequenced organism
* Reduce the huge number of small read sequences into larger, more manageable genomic contigs

Required inputs
---------------

* Sequenced reads in gzipped fastq format.

  * each sample is represented by two gzipped fastq files
  * standard output files of paired-end sequencing

::

   |-- reads/original
           |-- <sample_1>_R1.fastq.gz
           |-- <sample_1>_R2.fastq.gz
           |-- <sample_2>_R1.fastq.gz
           |-- <sample_2>_R2.fastq.gz

Generated outputs
-----------------

* Assembled genomic sequences (contigs) in fasta format

* Reports to assess quality of assembly

* Graph visualisation of assembly to visually assess its complexity

Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/mhv

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_assembly.yaml

Example configuration:

.. literalinclude:: ../../example/mhv/config_assembly.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of preprocess and mapping with the `MultiQC <https://multiqc.info/>`_
* Connect contigs into scaffolds based on known genomic sequence of related organism
* Aggregate quast results of individual samples into summary report


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report
   /pipelines/preprocess