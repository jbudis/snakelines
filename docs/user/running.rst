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
* Miniconda (tested on 4.5.11)


Installation
---------------

Sources codes for SnakeLines pipelines are stored at `GitHub repository <https://github.com/jbudis/snakelines>`_.
You may download them directly or clone them using git.

.. code:: bash

   # Download directly
   wget https://github.com/jbudis/snakelines/archive/master.zip
   unzip master.zip

   # Or clone using git
   git clone https://github.com/jbudis/snakelines

Compiling is not required, scripts are ready for use right after download.

Directory structure of input files
----------------------------------

SnakeLines takes as input files of two types; reads and reference-based files.
Read files in Fastq files represent sequenced DNA fragments.
SnakeLines in current version supports only paired-end reads that are represented by pair of files, with _R1.fastq.gz and _R2.fastq.gz extension.
Read files must be gzipped and located in the reads/original directory.

Reference-based files define sequence and annotation of reference genome.
Sequence file must be in the fasta format and have to be located in the reference/<genome>/<genome>.fa file.
Annotation is typically bed file with targeted genomic regions and have to be located in the reference/<genome>/<genome>/annotation/<panel>/regions.bed file.

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

Bioinformatic tools typically require preprocessed reference sequences to condensed files called indices.
All required reference indices and auxiliary files are generated, when necessary, during pipeline execution.

Reference directories with frequently used references may be linked to the project directory, to avoid redundant copies and repeated creation of sequence indices.
For example, if you have fasta file for human genome in separate directory (/data/genome/human/hg38-ucsc/hg38.fa), you may link it to example project (/data/projects/example) using

.. code:: bash

   ln --symbolic \
      /data/genome/human/hg38-ucsc \
      /data/projects/example/reference/hg38

Make sure, that the name of the link is the same as the name of the fasta file (without .fa suffix).

Running scripts
---------------

All SnakeLines pipelines are defined only by their configuration file in human-readable yaml format.
We recommend to copy the configuration file into the project directory.
This way, configuration for the pipeline is project specific, and therefore would not be shared between different projects.

Example project structure with configuration file copied from the <snakelines_dir>/example/mhv/
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
   |-- config_variant_calling.yaml

Edit config_variant_calling.yaml file according to your preference.
Each configured attribute is explained by a comment in the file.

Now you may run SnakeLines pipeline using Snakemake.
You need to specify one additional attribute, to tell Snakemake, where are SnakeLines sources located.
For example, if SnakeLines sources has been downloaded to the /usr/local/snakelines directory, use:

.. code:: bash

   snakemake \
      --snakefile /usr/local/snakelines/snakelines.snake \
      --configfile config_vairant_calling.yaml

Snakemake is very flexible in workflow execution, see `detailed documentation <https://snakemake.readthedocs.io/en/stable/executable.html#all-options>`_ and `useful bash aliases for SnakeLines <./aliases.html>`_.
