Analyse variants from nanopore reads
====================================

Apply the variant calling pipeline to reads from Nanopore sequencer. 

Purpose
-------

* Analyse variants in Nanopore reads

Required inputs
---------------

* Single-end reads from Nanopore sequencer in gzipped fastq format

  * Each sample is represented by a single gzipped fastq file

* Reference genome in fasta format

::

   |-- reads/original
           |-- <sample_1>.fastq.gz
           |-- <sample_2>.fastq.gz
   |-- reference/<reference>
           |-- <reference>.fa

Generated outputs
-----------------

* List of identified variants in VCF file, filtered by user-defined criteria

* Summary PDF report to assess quality of reads, mapping and variant calling


Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/nanopore

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config.yaml \
      --use-conda

Example configuration:

.. literalinclude:: ../../example/nanopore/config.yaml
   :language: yaml

