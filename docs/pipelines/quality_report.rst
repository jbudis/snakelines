Read quality report
===================

Generate summary HTMl report with quality statistics of sequenced reads, including

* Overall sequence quality
* Per base sequence quality
* GC content
* Sequence length
* Sequence duplication
* Adapter content

The pipeline uses `FastQC utility <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_ to generate quality reports for individual samples.
Individual reports are aggregated into the summary HTML report using custom scripts.

Purpose
-------

* Quickly assess quality of sequencing run
* Identify potential problems with downstream analysis - avoid sequencing artefacts
* Important to properly set configuration of downstream trimming analysis

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

* Summary HTML table with quality statistics of sequenced reads of multiple samples
* Individual FastQC reports

Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/mhv

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_quality_report.yaml

Example configuration:

.. literalinclude:: ../../example/mhv/config_quality_report.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of multiple samples with the `MultiQC <https://multiqc.info/>`_
