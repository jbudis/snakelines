Analyse methylation profiles
============================

Identify genomic regions with and without methylation.
The pipeline expects reads with the bisulfide conversion.

Purpose
-------

* Epigenetic marker for several diseases (e.g. oncology)
* Compare between samples with different phenotype (e.g. tissues)

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

* Summary report of methylation profiles in sequenced samples


Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/genomic

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_methylseq.yaml \
      --use-conda

Example configuration:

.. literalinclude:: ../../example/genomic/config_methylseq.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of preprocess and mapping with the `MultiQC <https://multiqc.info/>`_
* Include coverage tracks (Bismark can produce them as well)


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report
   /pipelines/preprocess
   /pipelines/mapping