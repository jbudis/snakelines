Analyse single-end Illumina samples
===================================

Apply the variant calling pipeline to single-end reads from Illumina sequencer. 

Purpose
-------

* Sometimes we need to process Illumina reads as single-end
* Single-end sequencing is cheaper than paired-end

Required inputs
---------------

* Single-end reads from Illumina sequencer in gzipped fastq format

  * Each sample is represented by a single gzipped fastq file
  * Half of the standard output files of paired-end sequencing

* Reference genome in fasta format

::

   |-- reads/original
           |-- <sample_1>_R1.fastq.gz
           |-- <sample_2>_R1.fastq.gz
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

   cd /usr/local/snakelines/example/illumina_sinle_end

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config.yaml \
      --use-conda

Example configuration:

.. literalinclude:: ../../example/illumina_single_end/config.yaml
   :language: yaml

