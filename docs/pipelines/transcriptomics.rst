Analyse gene expression
=======================

Determine level of expression of individual transcripts from sequenced RNA.
Compare expressions between samples from two conditions to identify transcripts with changed expression.

Purpose
-------

* Identify which transcripts are expressed
* Compare expression of transcripts with each other in a single sample
* Compare changes of expression between samples with different conditions, e.g.

  * environment
  * tissues
  * phenotypic trait

* Assess impact of a gene on biochemical pathways (compare wild type with sample without the gene)

Required inputs
---------------

* Sequenced reads in gzipped fastq format.

  * each sample is represented by two gzipped fastq files
  * standard output files of paired-end sequencing

* Transcripts of reference genome in fasta format

* Sample metadata file in TSV format

  * each row describe a single sequenced sample
  * the first, 'sample-id' column represents names of fastq files without _R[12].fastq.gz part, e.g. sample_1, sample_2
  * other columns represent attributes of sequenced samples that may be used to separate samples into categories for comparison, e.g. environment, tissue

::

   |-- description
           |-- sample-metadata.csv
   |-- reads/original
           |-- <sample_1>_R1.fastq.gz
           |-- <sample_1>_R2.fastq.gz
           |-- <sample_2>_R1.fastq.gz
           |-- <sample_2>_R2.fastq.gz
   |-- reference/<reference>
           |-- <reference>.transcripts.fa

Generated outputs
-----------------

* Summary table of number of transcripts per sample, normalised by read count and length of transcripts

* Graphical, 2D PCA comparison of samples based on their expression to visually assess relationships between them

* Set of transcripts with significant set of expression between selected conditions


Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/genomic

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_transcriptomics.yaml

Example configuration:

.. literalinclude:: ../../example/rnaseq/config_transcriptomics.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of preprocess and mapping with the `MultiQC <https://multiqc.info/>`_
* Add mapped based approach with resulting coverage tracks


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report
   /pipelines/preprocess
