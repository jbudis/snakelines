How to write documentation
==========================

SnakeLines documentation is written in the `rst format <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ (ReStructuredText).
Documentation pages are generated from

* static .rst files located in the "docs/development/" and "docs/user/" directories
* automatically generated .rst files according to docstrings of individual rules and pipelines (directories "docs/rules/" and "docs/pipeline/")

Document each new rule
----------------------

The documentation for each rule is generated from the docstring of each rule. The docstring follows the ReStructuredText (reST) Docstring style from python (`Examples of different docstring styles for python <http://queirozf.com/entries/docstrings-by-example-documenting-python-code-the-right-way>`_), but recognizes only three basic directives:

.. code-block:: bash

   :input <input_name>: <input_description>
   :output <output_name>: <output_description>
   :param <param_name>: <parameter_description>

An example of a correct documentation of a rule:

.. code-block:: bash

   rule seqtk__subsample_reads:
       """
       Randomly select user-configured number of reads from fastq files.
       :input r1: filename with R1 reads
       :input r2: filename with R2 reads
       :output r1: filename with subsampled R1 reads
       :output r2: filename with subsampled R2 reads
       :param seed: int - seed of the random generator for subsampling
       :param n_reads: int - number of reads to keep in subsampled set
       """

The generated documentation corresponding to this example is in:

.. toctree::
   rule_example.rst


Build documentation locally
---------------------------

To regenerate all automatically generated .rst files, you should run the following command in the root directory of SnakeLines (you need to have python installed on your system):

.. code-block:: bash

   python src/documentation/documentation.py

You may need to install `Sphinx tool <http://www.sphinx-doc.org/en/master/>`_ first.

.. code-block:: bash

   pip install sphinx

Use prepared Makefile to convert .rst files to .html pages (Windows users need to have `"make" for Windows <http://gnuwin32.sourceforge.net/packages/make.htm>`_ installed)

.. code-block:: bash

   cd docs/
   make html

Html pages are stored in the _build/html directory.
You may review them using your favorite browser, for example

.. code-block:: bash

   firefox _build/html/index.html


Publish documentation
---------------------

To apply local changes in the documentation to the public `ReadTheDocs web documentation <https://snakelines.readthedocs.io/en/latest/>`_ (you are probably reading it now), just push docs/ directory to the github master branch.
At first, the ReadTheDocs server will be notified by the Github server, that changes were made.
Then, the master branch of the SnakeLines repository would be cloned to the ReadTheDocs server and built using Sphinx server.
Finally, generated Html files would be published.

The process takes a while, so you will see changes after few minutes.


Useful links
------------

* `How to set-up documentation from scratch <https://dont-be-afraid-to-commit.readthedocs.io/en/latest/documentation.html>`_ - more comprehensive version of this page