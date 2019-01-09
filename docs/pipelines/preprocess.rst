Preprocess sequencing reads
===========================

Remove parts or whole reads that are artifacts of laboratory sequencing process.
They may blur a downstream analysis, and so lead to incorrect conclusions in their interpretations.
According to provided configuration, preprocess may include:

* Removal of low quality parts of reads or adapters
* Removal of PCR duplicates
* Select a fixed number of reads from each sample to ensure consistency
* Removal of contamination from a known genome
* Selecting only fragments from a known genome
* Merging paired reads based on their read sequence overlap

Purpose
-------

* Remove sequencing artifacts
* Clean-up sequencing data for downstream analysis
* Avoid false interpretation of data analysis results due to laboratory-induced or technical bias


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

* Refined reads in gzipped fastq format.

  * each sample is represented by two gzipped fastq files

* Quality reports for resulting and intermediate reads to assess effect of individual preprocess steps

Example
-------

How to run example:

.. code-block:: bash

   cd /usr/local/snakelines/example/mhv

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_preprocess.yaml

Example configuration:

.. literalinclude:: ../../example/mhv/config_preprocess.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of multiple samples and processing steps with the `MultiQC <https://multiqc.info/>`_


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report