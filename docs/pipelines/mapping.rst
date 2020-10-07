Map reads to a reference genome
===============================

Align sequenced reads to a reference genome to find the most probable genomic positions of origin.
Reference genome may represent DNA of a single organism, multiple organisms, transcripts...

Purpose
-------

* Identify the source of sequenced material
* Assess composition of sequenced material
* Removal of a known contamination
* Purification of sequenced reads - select only reads from the target organism
* Essential step for several downstream analysis types - variant calling, transcriptomics ...


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

* Mapped reads in sorted, indexed BAM format with marked duplicates

* Reports to assess mapping quality

  * individual report for each sample
  * summary report for comparison of multiple samples

Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/genomic

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_mapping.yaml

Example configuration:

.. literalinclude:: ../../example/genomic/config_mapping.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of preprocess and mapping with the `MultiQC <https://multiqc.info/>`_
* Realignment postprocess step to refine alignment in indel regions


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report
   /pipelines/preprocess