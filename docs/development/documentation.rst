How to write documentation
==========================

SnakeLines documentation is written in the `rst format <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ (ReStructuredText).
Documentation pages are generated from

* static .rst files located in the docs/ directory
* automatically generated .rst files according to docstrings of individual rules and pipelines.

Build documentation locally
---------------------------

You may need to install `Sphinx tool <http://www.sphinx-doc.org/en/master/>`_ tool first.

.. code-block:: bash

   pip install sphinx

Use prepared Makefile to convert .rst files to .html pages

.. code-block:: bash

   cd docs/
   make html

Html pages are stored in the _build/html directory.
You may review them using your favorite browser, for example

.. code-block:: bash

   firefox _build/html/index.html


Publish documentation
---------------------

To apply local changes in the documentation to the public `ReadTheDocs web documentation <https://snakelines.readthedocs.io/en/latest/>`_ (you are probably reading it not), just push docs/ directory to the github master branch.
At first, the ReadTheDocs server will be notified by the Github server, that changes were made.
Then, the master branch of the SnakeLines repository would be cloned to the ReadTheDocs server and built using Sphinx server.
Finally, generated Html files would be published.

The process takes a while, so you will see changes after few minutes.


Useful links
------------

* `How to set-up documentation from scratch <https://dont-be-afraid-to-commit.readthedocs.io/en/latest/documentation.html>`_ - more comprehensive version of this page