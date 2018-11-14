How to write documentation
==========================

SnakeLines documentation is written in the `rst format <http://docutils.sourceforge.net/docs/user/rst/quickstart.html>`_ (ReStructuredText).
Documentation pages are generated from

* static .rst files located in the docs/ directory
* automatically generated .rst files according to docstrings of individual rules and pipelines.

Html pages are generated using `Sphinx tool <http://www.sphinx-doc.org/en/master/>`_.

Build documentation locally
---------------------------

You may need to install sphinx tool first.

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

Useful links
------------

* `How to set-up documentation from scratch <https://dont-be-afraid-to-commit.readthedocs.io/en/latest/documentation.html>`_ - more comprehensive version of this page