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

Python3 module dependencies
---------------------------

Snakelines requires several Python3 modules to be installed prior to Snakelines. All of these modules can be installed via Conda package manager or via Pip3 with the exception of Tkinter that can be obtained via e.g. apt-get.

Note that Snakelines requires specific versions of these modules.
List of required Python3 modules:

* oyaml=0.7
* pandas=0.19.2
* biopython=1.72
* seaborn=0.9.0
* bs4=0.0.1
* weasyprint=0.3
* pysam=0.15.1
* openpyxl=2.5.12
* scikit-learn=0.18
* tk=3.6.7-1

Snakelines will be later shipped as a Conda package in Anaconda repository.

Conda channel dependencies
--------------------------

In order for Snakelines Conda virtual enviroments to work, user has to add several Anaconda repository channels to Anaconda using these commands:

.. code:: bash

   conda config --add channels bioconda
   conda config --add channels g2554711
   conda config --add channels g2554711/label/bioconda
   conda config --add channels conda-forge
   conda config --add channels agbiome
   conda config --add channels rsmulktis
   conda config --add channels moustik
   
``--use-conda`` option in ``Snakemake`` command will enable use of predefined virtual enviroments in Snakelines.
   
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
      --configfile config_variant_calling.yaml

Snakemake is very flexible in workflow execution, see `detailed documentation <https://snakemake.readthedocs.io/en/stable/executable.html#all-options>`_ and `useful bash aliases for SnakeLines <./aliases.html>`_.

Multi-threading
---------------

SnakeLines executes tools that support parallelization on multiple cores, using standard `Snakemake features <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads>`_.
The number of threads for each task may be specified in a Snakemake call as:

.. code:: bash

   snakemake \
      --snakefile /usr/local/snakelines/snakelines.snake \
      --configfile config_variant_calling.yaml \
      --config threads=8

Alternately, user may specify the number of threads directly in a configuration file:

.. code:: yaml

   threads: 16                         # Number of threads to use in analysis
   samples:                            # List of sample categories to be analysed
       - name: example.*               # Regex expression of sample names to be analysed (reads/original/example.*_R1.fastq.gz)
         reference: mhv                # Reference genome for reads in the category (reference/mhv/mhv.fa)

   report_dir: report/public/01-assembly   # Generated reports and essential output files would be stored there

   reads:                              # Prepare reads and quality reports for downstream analysis
       preprocess:                     # Pre-process of reads, eliminate sequencing artifacts, contamination ...

           trimmed:                    # Remove low quality parts of reads
               method: trimmomatic     # Supported values: trimmomatic
               temporary: False        # If True, generated files would be removed after successful analysis
               crop: 500               # Maximal number of bases in read to keep. Longer reads would be truncated.
               quality: 20             # Minimal average quality of read bases to keep (inside sliding window of length 5)
               headcrop: 20            # Number of bases to remove from the start of read
               minlen: 35              # Minimal length of trimmed read. Shorter reads would be removed.

SnakeLines uses 1 core by default, if the number of threads is not specified.