Analyse microbial composition
=============================

Determine composition of organisms (typically microbial) in sequenced samples.
Data analysis pipeline expect data from shotgun sequencing of a single target gene, for example 16S.
Visually assess changes in compositions between several samples.

Purpose
-------

* Identify which organisms are present in sequenced samples
* Compare an abundance of organisms with each other in a single sample
* Compare abundances of organisms between samples

Required inputs
---------------

* Sequenced paired-end reads from Illumina sequencer in gzipped fastq format.

  * each sample is represented by two gzipped fastq files
  * standard output files of paired-end sequencing

* .tax file with taxonomic description of each organism in reference fasta

  * TSV format
  * no header
  * each row represents a single target gene from the reference fasta
  * the first column represents the identifier of the gene (the first word after > symbol)
  * the second column represents taxonomy of the organism
  * for example

::

   KF494428.1.1396 Bacteria;Epsilonbacteraeota;Campylobacteria;Campylobacterales;Thiovulaceae;Sulfuricurvum;Sulfuricurvum sp. EW1;;;;;;;;
   AF506248.1.1375 Bacteria;Cyanobacteria;Oxyphotobacteria;Nostocales;Nostocaceae;Nostoc PCC-73102;Nostoc sp. 'Nephroma expallidum cyanobiont 23';;;;;;;;
   EF603722.1.1487 Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Muribaculaceae;;;;;;;;;;
   ...

::

   |-- reads/original
           |-- <sample_1>_R1.fastq.gz
           |-- <sample_1>_R2.fastq.gz
           |-- <sample_2>_R1.fastq.gz
           |-- <sample_2>_R2.fastq.gz


In case, you are using Metaxa2 classifier, you also need reference genomes of target genes in fasta format.

::

   |-- reference/<reference>
           |-- <reference>.fa
           |-- <reference>.tax

Generated outputs
-----------------

* Summary table with the number of reads per organism in sequenced sample

* Pie plot with hierarchical taxonomical composition of the samples

* Barplots for visual comparison of organisms abundance across sequenced samples


Example
-------

Metaxa2 classifier
~~~~~~~~~~~~~~~~~~

How to run example:

You need to provide reference sequences for this example to run.
Sequences of the 16S gene (store them as reference/silva-16S/silva-16S.fa) may be downloaded from the `Silva database <https://www.arb-silva.de/>`_.
Sequences of the ITS gene (store them as reference/unite/unite.fa) may be downloaded from the `Unite database <https://unite.ut.ee//>`_.

.. code-block:: bash

   cd /usr/local/snakelines/example/metagenomics

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_metagenomics.yaml \
      --use-conda

Example configuration:

.. literalinclude:: ../../example/metagenomic/config_metagenomic.yaml
   :language: yaml

RDP classifier
~~~~~~~~~~~~~~

How to run example:

RDP classifier has already prebuilt databases that may be used, specifically 16srrna, fungallsu, fungalits_warcup, fungalits_unite


.. code-block:: bash

   cd /usr/local/snakelines/example/metagenomics

   snakemake \
      --snakefile ../../snakelines.snake \
      --configfile config_metagenomics_rdp.yaml

Example configuration:

.. literalinclude:: ../../example/metagenomic/config_metagenomic_rdp.yaml
   :language: yaml

Planned improvements
--------------------

* Aggregate quality statistics of preprocess and mapping with the `MultiQC <https://multiqc.info/>`_
* Add toy reference genome to example
* Qiime2 based analysis


Included pipelines
------------------

.. toctree::
   :maxdepth: 2

   /pipelines/quality_report
   /pipelines/preprocess_paired_end
