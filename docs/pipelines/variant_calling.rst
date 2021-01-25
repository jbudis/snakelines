Call variants
=============

Identify alternations of a genomic material of a sequenced individual with respect to a model, reference genome.
That includes small point mutations (SNPs), small insertions and deletions.

Purpose
-------

* Identify genomic mutation that causes phenotypic trait of interest
* Identify genomic differences between close species
* Important for various genetic tests

  * inherited diseases
  * de novo mutations
  * oncological diseases - screening and assessing its type

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

* List of identified variants in VCF file, filtered by user-defined criteria

* Summary PDF report to assess quality of reads, mapping and variant calling


Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/genomic

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_variant_calling.yaml \
      --use-conda

Example configuration:

.. literalinclude:: ../../example/genomic/config_variant_calling.yaml
   :language: yaml

Planned improvements
--------------------

* Call variants only in pre-defined regions (BED file)
* Aggregate quality statistics of preprocess and mapping with the `MultiQC <https://multiqc.info/>`_
* Annotation of variants
* Filtering of variants based on external annotations
* Interpretation of variants


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report
   /pipelines/preprocess
   /pipelines/mapping