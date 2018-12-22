Source code structure
=====================

`Downloaded source codes <../user/running.html#installation>`_ are organised in following directories:

* docs - Documentation for the SnakeLines project (you are reading it now)
* example - Small toy examples of read files and references to try SnakeLines pipelines
* legacy - Old source codes from the first version of SnakeLines, waiting to be included
* **rules** - Hierarchy of simple Snakemake rules
* src - SnakeLines specific source codes that facilitate particular nature of its configuration and source code structure

Users should modify only config yaml files, which define "what would be generated".
Developers would be more interested in the rules directory, with logic behind the files creation; "how it would be generated".