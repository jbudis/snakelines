Download reference and map reads
================================

At first, download sequences from NCBI according to selected genbank ids and prepare reference database.
Then proceed with the standard mapping pipeline.
Align sequenced reads to a reference genome to find the most probable genomic positions of origin.
Reference genome may represent DNA of a single organism, multiple organisms, transcripts...

The reference part of the configuration may be similarly used for other types of analysis, such as variant calling, methylation...


Purpose
-------

* Download genomic sequences and prepare reference database
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

::

   |-- reads/original
           |-- <sample_1>_R1.fastq.gz
           |-- <sample_1>_R2.fastq.gz
           |-- <sample_2>_R1.fastq.gz
           |-- <sample_2>_R2.fastq.gz

Generated outputs
-----------------

* Reference sequence database with taxonomies

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
      --configfile config_download_reference_and_mapping.yaml \
      -- use-conda

Example configuration:

.. literalinclude:: ../../example/genomic/config_download_reference_and_mapping.yaml
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
   /pipelines/mapping