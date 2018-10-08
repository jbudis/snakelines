Running SnakeLines
==================

SnakeLines are designed to simplify running of complex bioinformatics pipelines without manual install of any additional libraries and bioinformatics tools.
All dependencies are installed using Conda into a new virtual environment.
This way, there should be no conflicts between already installed tools and SnakeLines specific configuration.



Software requirements
---------------------

Minimal software requirements are:

* Linux (tested on Ubuntu 16.04.1)
* `SnakeMake <https://snakemake.readthedocs.io/en/stable/>`_ (tested on 5.2.2)

Directory structure of input files
----------------------------------

SnakeLines takes as input files of two types.
Read files in Fastq files represent sequenced DNA fragments.
SnakeLines in current version supports only paired-end reads that are represented by pair of files, with _R1.fastq.gz and _R2.fastq.gz extension.
Read files must be gzipped and located in the reads/original directory.

Reference-based files define sequence and annotation of reference genome.
Sequence file must be in the fasta format and have to be located in the reference/<genome>/<genome>.fa file.
Annotation is typically bed file with targeted genomic regions and habe to be located in the reference/<genome>/<genome>/annotation/<panel>/regions.bed file.

Example minimal input file configuration for a reference (hg38), targeted panel (sureselect6) and reads for samples example_A and example_B:
::
   |-- reads/original
           |-- example_A_R1.fastq.gz
           |-- example_A_R2.fastq.gz
           |-- example_B_R1.fastq.gz
           |-- example_B_R2.fastq.gz
   |-- reference/hg38
           |-- hg38.fa
           |-- annotation/sureselect6
                   |-- regions.bed

