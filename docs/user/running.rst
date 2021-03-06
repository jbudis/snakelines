Running SnakeLines
==================

SnakeLines are designed to simplify running of complex bioinformatics pipelines without manual install of any additional libraries and bioinformatics tools.
All dependencies are installed using Conda into a new virtual environment.
This way, there should be no conflicts between already installed tools and SnakeLines specific configuration.



Software requirements
---------------------

Minimal software requirements are:

* Linux (tested on Ubuntu 20.04)
* Miniconda (tested on 4.8.3)
   
Installation
---------------

There are two ways for SnakeLines to be installed: as Conda package and by cloning the Github repository. We recommend the first option.

Installation as a Conda package
--------------------------------

.. code:: bash
   
   # In order for SnakeLines Conda virtual enviroments to work, user has to add several Anaconda repository channels to Conda.
   
   conda config --add channels bioconda
   conda config --add channels g2554711
   conda config --add channels g2554711/label/bioconda
   conda config --add channels conda-forge
   conda config --add channels agbiome
   conda config --add channels rsmulktis
   conda config --add channels moustik
  
   conda create --name snakelines-env -c bioconda snakelines
   conda activate snakelines-env
   

Installation from Github repository
------------------------------------

Snakelines requires several Python3 modules to be installed prior to Snakelines. All of these modules can be installed via Conda package manager or via Pip3 with the exception of Tkinter that can be obtained via e.g. apt-get.

.. code:: bash

    pip install numpy==1.19.2 oyaml==0.9 pandas==1.1.3 biopython==1.78 seaborn==0.11.0 bs4==4.9.3 weasyprint==51 pysam==0.16.0.1 openpyxl==3.0.5 scikit-bio==0.5.6 jinja2==2.11.2 snakemake==5.13.0

Sources codes for SnakeLines pipelines are stored at `GitHub repository <https://github.com/jbudis/snakelines>`_.
You may download them directly or clone them using git.

.. code:: bash

   # Download directly
   wget https://github.com/jbudis/snakelines/archive/master.zip
   unzip master.zip

   # Or clone using git
   git clone https://github.com/jbudis/snakelines

Compiling is not required, scripts are ready for use right after download.

Note that GATK3 has to be installed separately (if the user intends to use it) as the license prohibits GATK3 to be distributed via Conda repository.

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

Running scripts
---------------

All SnakeLines pipelines are defined only by their configuration file in human-readable yaml format.
We recommend to copy the configuration file into the project directory.
This way, configuration for the pipeline is project specific, and therefore would not be shared between different projects.

Example project structure with configuration file copied from the <snakelines_dir>/example/genomic/
::

   |-- reads/original
           |-- example_A_R1.fastq.gz
           |-- example_A_R2.fastq.gz
           |-- example_B_R1.fastq.gz
           |-- example_B_R2.fastq.gz
   |-- reference/mhv
           |-- mhv.fa
   |-- config_variant_calling.yaml

Edit config_variant_calling.yaml file according to your preference.
Each configured attribute is explained by a comment in the file.

Now you may run SnakeLines pipeline using Snakemake.
You need to specify one additional attribute, to tell Snakemake, where are SnakeLines sources located.
If SnakeLines was installed as a Conda package, there is a wrapper script available, so the resulting command will be:

.. code:: bash

   snakelines --configfile config_variant_calling.yaml --use-conda

You can set all parameters from ``snakemake`` in ``snakelines`` wrapper.

In case of installation from Github repository, if SnakeLines sources have been downloaded to the /usr/local/snakelines directory, use:

.. code:: bash

   snakemake \
      --snakefile /usr/local/snakelines/snakelines.snake \
      --configfile config_variant_calling.yaml \
      --use-conda

Snakemake is very flexible in workflow execution, see `detailed documentation <https://snakemake.readthedocs.io/en/stable/executable.html#all-options>`_ and `useful bash aliases for SnakeLines <./aliases.html>`_.


Reference files
---------------

Bioinformatic tools typically require preprocessed reference sequences to condensed files called indices.
All required reference indices and auxiliary files are generated, when necessary, during pipeline execution.

Reference directories with frequently used references may be linked to the project directory, to avoid redundant copies and repeated creation of sequence indices.
For example, if you have fasta file for human genome in separate directory (/data/genome/human/hg38-ucsc/hg38.fa), you may link it to example project (/data/projects/example) using

.. code:: bash

   ln --symbolic \
      /data/genome/human/hg38-ucsc \
      /data/projects/example/reference/hg38

Make sure, that the name of the link is the same as the name of the fasta file (without .fa suffix).

Download sequences from NCBI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SnakeMake can prepare reference database from the provided list of genbank ids.
At first, sequences with specified ids in the configuration file would be downloaded from NCBI and aggregated into a single fasta file.
Next, .tax file with taxonomies of downloaded sequences will be created.
Finally, created sequence and taxonomy files would be used as the reference for downstream analysis.

See example pipeline for `the mapping with downloaded reference <../pipelines/download_reference_and_mapping.html>`_.
Other pipelines may be updated accordingly, you just need to include the ``reference`` block of configuration:

.. code:: yaml

   reference:
      download:
         method: entrez               # Supported values: entrez
         email: FILLME@SOMEMAIL.COM   # Inform NCBI who you are to contact you in case of excessive use.
         mhv_ncbi:                    # List of genbank ids to download, one list for each reference database
            - U97553.2
            - AF127083.1



Use reference indices without fasta
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, it makes sense to keep only indices, without the primary fasta file.
For example, huge genomic databases provided by NCBI are already packed into Blast indices.
Downloading original fasta files and generating indices is a huge burden to memory and computational capacity of a cluster.

In such scenarios, you may use just downloaded indices, without the primary fasta file.
Keep in mind, that such reference could be used only for one tool, as Blast in this example.
Assuming, you downloaded Blast indices are stored at
::

   |-- /data/genome/metagenome/blast/nt/
           |-- nt.00.nhd
           |-- nt.00.nhi
           |-- nt.00.nhr
           |-- nt.01.nhd
           |-- nt.01.nhi
           |-- nt.01.nhr
           |-- ...
           |-- nt.60.nhd
           |-- nt.60.nhi
           |-- nt.60.nhr
           |-- nt.nal
           |-- taxdb.btd
           |-- taxdb.bti


You may link index directly to the project using

.. code:: bash

   ln --symbolic \
      /data/genome/metagenome/blast/nt/ \
      /data/projects/example/reference/nt/blast_index

Such databases should be labelled with ``prebuilt: True`` value in the configuration, to avoid validation messages for missing fasta file:

.. code:: bash

   samples:                           # List of sample categories to be analysed
      - name: .*-16S                  # Regex expression of sample names to be analysed (reads/original/.*-16S_R1.fastq.gz)
        reference: 16srrna            # RDP classifier Supported values: 16srrna, fungallsu, fungalits_unite, fungalits_warcup
        prebuilt: True                # Reference sequence reference/{reference}/{reference}.fa does not exist, but all required indices are already prepared

Multi-threading
---------------

SnakeLines executes tools that support parallelization on multiple cores, using standard `Snakemake features <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads>`_.
The number of threads for each task may be specified in a Snakemake call as:

.. code:: bash

   snakemake \
      --snakefile /usr/local/snakelines/snakelines.snake \
      --configfile config_variant_calling.yaml \
      --use-conda \
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
