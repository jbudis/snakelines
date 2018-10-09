
Why use SnakeLines
====================

Several workflow management systems have been proposed to date, most notably `Galaxy <https://galaxyproject.org/>`_ and `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_.
SnakeLines extends traditional SnakeMake rules with following benefits:


Ready-to-use pipelines
~~~~~~~~~~~~~~~~~~~~~~

SnakeLines contains wide set of computational pipelines that are ready to use immediately after downloading SnakeLines sources.
In addition to standard secondary analysis, such as mapping and assembly, SnakeMake facilitate:

* identification of variation in mapped reads
* identification of viruses in reads or assembled contigs
* classify assembled contigs
* metagenomics analysis
* transcriptomics (in preparation)
* chromatin-seq (in preparation)


Reporting
~~~~~~~~~

Each step of analysis is supported by graphical or tabular reports.
Essential output and report files are stored in separate, report directory.
In addition to reports that are generated for each sample individually, SnakeLines store also aggregated reports to examine quality of several samples at once.


No need to install additional tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SnakeLines requires only Snakemake and Miniconda to be installed.
The rest of bioinformatic tools required for analysis are compiled automatically.
Tools are installed into separated virtual environments.
This way, they do not break any installation or dependencies on the operating system.


Easily configurable
~~~~~~~~~~~~~~~~~~~

Each pipeline is parametrized in a single, YAML based configuration file.
Configuration is useful overview over analysis, since every step of the analysis is configured there in sequential order.
You may easily swap tools or change their parameters by adjusting the configuration file.


Easily extensible
~~~~~~~~~~~~~~~~~

Structure of the SnakeLines sources allows to easily extend existing pipelines, or just replace method in some step with custom solution.
