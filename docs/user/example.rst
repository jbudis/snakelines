Quick start
===========

This guide will show you the basic usage of SnakeLines pipelines on small toy example.
We assume that you have Linux system installed (tested on Ubuntu 20.04).

Install required dependencies
-----------------------------

Follow installation guides for `Miniconda <https://conda.io/docs/user-guide/install/index.html>`_.
There is no need to install any other tools manually.
All required bioinformatic tools are installed automatically by SnakeLines inside dedicated virtual environments to avoid dependency conflicts with already installed tools.

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
   

Alternately you may `download sources directly <running.html#installation>`_.

Execute pipeline
----------------

Toy example read files and references are stored at the `/example` directory of downloaded sources.
Assuming, that sources are downloaded in /usr/local/snakelines, `variant_calling` pipeline may be executed using

.. code-block:: bash

   conda activate snakelines-env
   
   cd /usr/local/snakelines/example/genomic

   snakemake \
      --configfile config_variant_calling.yaml \
      --use-conda
