For Developers
==============

Installation
------------

To use Lumos, first install it in editable mode using pip:

.. code-block:: console

   (.venv) $ pip install lumos-sat -e

Building Lumos
--------------

To build Lumos, use the following commands:

.. code-block:: console

   (.venv) $ python build
   (.venv) $ twine check dist/*
   (.venv) $ twine upload -r testpypi dist/*
   (.venv) $ twine upload dist/*

Building Documentation Locally
------------------------------

.. code-block:: console

   (.venv) $ cd docs
   (.venv) $ make clean
   (.venv) $ make html
   (.venv) $ start build/html/index.html