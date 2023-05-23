Installation from Source
------------------------

First, clone the Github:

.. code-block:: console

   $ git clone https://github.com/Forrest-Fankhauser/lumos-sat.git

Next, create a virtual environment with the necessary packages:

.. code-block:: console

   $ python -m venv venv
   $ venv/Scripts/activate.bat
   (venv) $ pip install -r requirements.txt
   (venv) $ pip install -r dev_requirements.txt

Finally, install Lumos in editable mode:

.. code-block:: console

   $ cd lumos
   $ python -m pip install -e .
   
Building a Lumos Distribution
-----------------------------

First update the Lumos version number in pyproject.toml

To build Lumos, use the following commands:

.. code-block:: console

   (.venv) $ python -m build
   (.venv) $ twine check dist/*

First check your distribution thoroughly using Test PyPi.

.. code-block:: console

   (.venv) $ twine upload -r testpypi dist/*

If that works, you can upload to PyPi

.. code-block:: console

   (.venv) $ twine upload dist/*

Building Documentation Locally
------------------------------

.. code-block:: console

   (.venv) $ cd docs
   (.venv) $ make clean
   (.venv) $ make html
   (.venv) $ start build/html/index.html