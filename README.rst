SCAKS
=====

SCAKS: **S**\ hanghai **CA**\ talysis **K**\ inetics **S**\ oftware

Install
-------

1. Via pip (Not uploaded to PyPI yet):

   .. code:: shell

      pip install scaks

2. From source

   .. code:: shell

      git clone https://github.com/PytLab/scaks.git

      cd scaks
      python setup.py install

Test
----

.. code:: shell

   python setup.py test

or

.. code:: shell

   python -m scaks.tests.test_all

Documentation
-------------

For more details, see SCAKS documentation:
https://scaks-docs.readthedocs.io/en/latest/

Code structure
--------------

::

   ├── LICENSE                # License file
   ├── MANIFEST.in            # Define the list of files to include in the distribution
   ├── README.rst             # This file
   ├── docs                   # Documentation directory
   │   ├── LICENSE
   │   ├── Makefile           # Makefile to build documentation
   │   ├── README.rst
   │   ├── build              # Built HTML documentation
   │   └── source             # Documentation source code
   ├── examples               # Usage examples
   ├── requirements.txt       # Python dependencies information
   ├── scaks                  # SCAKS main package
   │   ├── __init__.py        # Contains base class for microkinetic model components
   │   ├── compatutil.py      # A module for compatibility of Python2.x & Python3.x
   │   ├── correctors         # Energy correctors of microkinetic model
   │   ├── database           # Constant data information
   │   ├── descriptors        # Some descriptors definition used by model
   │   ├── errors             # Custom Python exception classes
   │   ├── functions.py       # Some common utility functions
   │   ├── models             # Kinetic models
   │   ├── mpicommons.py      # A wrapped MPI for Python interfaces used in SCAKS
   │   ├── parsers            # Parers of microkinetic model
   │   ├── plotters           # Plotters for energy profile plotting based on CatPlot
   │   ├── plugins            # On-the-fly analysis and hybrid method plugins and metaclasses
   │   ├── solvers            # Solvers used by microkinetic models
   │   ├── tests              # Unit test for SCAKS interfaces
   │   └── utilities          # Some utilities functions
   ├── scripts                # Some scripts for reference of users
   │   ├── mkm_plot           # Some plotting scripts
   │   └── mkm_run            # Run script
   ├── setup.cfg              # Python setup configure file
   └── setup.py               # Python setup script based on setuptools
