Variant detection in SARS-CoV-2
===============================

Process paired-end reads from the Illumina COVIDSeq protocol to characterize sequenced SARS-CoV-2 virus.

Purpose
-------

* detect genomic variants
* assign the sequenced SARS-CoV-2 into COVID clade
* make consensus sequence of the sequenced virus


Required inputs
---------------
* Sequences of primers used in PCR amplicon based enrichment (fasta format, artic.fa)
* Sequenced paired-end reads from Illumina sequencer in gzipped fastq format.

  * each sample is represented by two gzipped fastq files
  * standard output files of paired-end sequencing

* Reference genome in fasta format

::

   |-- description/adapters
           |-- artic.fa
   |-- reads/original
           |-- <sample_1>_R1.fastq.gz
           |-- <sample_1>_R2.fastq.gz
           |-- <sample_2>_R1.fastq.gz
           |-- <sample_2>_R2.fastq.gz
   |-- reference/<reference>
           |-- <reference>.fa


Generated outputs
-----------------

* Consensus sequence of the sequenced SARS-CoV-2 virus
* Detected genomic variants
* Clade assignment


Example
-------

How to run example:

Copy Fastq files from the Illumina COVIDSeq sequencing into /usr/local/snakelines/example/covidseq/reads/original.
Real-time examples can be also freely downloaded from the `ENA portal <https://www.ebi.ac.uk/ena/browser/view/PRJEB43444>`_.

.. code-block:: bash

   cd /usr/local/snakelines/example/covidseq

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_covidseq.yaml \
      --use-conda

Example configuration:

.. literalinclude:: ../../example/covidseq/config_covidseq.yaml
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
   /pipelines/preprocess_paired_end
