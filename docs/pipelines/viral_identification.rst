Identify viruses
================

Identify reads of viral origin in metagenomic samples.

Purpose
-------

* Identify specific strain of pathogenic virus in infected sample
* Assess virome composition of a sample

Required inputs
---------------

* Sequenced reads in gzipped fastq format.

  * each sample is represented by two gzipped fastq files
  * standard output files of paired-end sequencing

* Reference genome in fasta format

::

   |-- reads/original
           |-- <sample_1>_R1.fastq.gz
           |-- <sample_1>_R2.fastq.gz
           |-- <sample_2>_R1.fastq.gz
           |-- <sample_2>_R2.fastq.gz
   |-- reference/<reference>
           |-- <reference>.fa

Generated outputs
-----------------

* Hierarchical interactive pie plot (Krona) for visual assessment of viral composition in the sample


Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/mhv

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_viral_identification.yaml

Example configuration:

.. literalinclude:: ../../example/genomic/config_viral_identification.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of preprocess and mapping with the `MultiQC <https://multiqc.info/>`_
* Assembly based method


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report
   /pipelines/preprocess