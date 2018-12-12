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

   ├── LICENSE                           # License file
   ├── MANIFEST.in                       # Define the list of files to include in the distribution
   ├── README.rst                        # This file
   ├── docs                              # Documentation directory
   │   ├── LICENSE
   │   ├── Makefile                      # Makefile to build documentation
   │   ├── README.rst
   │   ├── build                         # Built HTML documentation
   │   └── source                        # Documentation source code
   ├── examples                          # Usage examples
   ├── requirements.txt                  # Python dependencies information
   ├── scaks                             # SCAKS main package
   │   ├── __init__.py                   # Contains base class for microkinetic model components
   │   ├── compatutil.py                 # A module for compatibility of Python2.x & Python3.x
   │   ├── correctors                    # Energy correctors of microkinetic model
   │   ├── database                      # Constant data information
   │   │   ├── __init__.py
   │   │   ├── elements_data.py
   │   │   ├── lattice_data.py
   │   │   └── thermo_data.py
   │   ├── descriptors                   # Some descriptors definition used by model
   │   │   ├── __init__.py
   │   │   ├── component_descriptors.py
   │   │   └── descriptors.py
   │   ├── errors                        # Custom Python exception classes
   │   ├── functions.py                  # Some common utility functions
   │   ├── models                        # Kinetic models
   │   │   ├── __init__.py
   │   │   ├── kinetic_model.py
   │   │   └── micro_kinetic_model.py
   │   ├── mpicommons.py                 # A wrapped MPI for Python interfaces used in SCAKS
   │   ├── parsers                       # Parers of microkinetic model
   │   │   ├── __init__.py
   │   │   ├── absolute_energy_parser.py
   │   │   ├── parser_base.py
   │   │   ├── relative_energy_parser.py
   │   │   └── rxn_parser.py
   │   ├── plotters                      # Plotters for energy profile plotting based on CatPlot
   │   │   ├── __init__.py
   │   │   ├── energy_profile_plotter.py
   │   │   └── plotter_base.py
   │   ├── plugins                       # On-the-fly analysis and hybrid method plugins and metaclasses
   │   │   ├── __init__.py
   │   │   ├── analysis.py
   │   │   ├── hybrid_methods.py
   │   │   └── metaclasses.py
   │   ├── solvers                       # Solvers used by microkinetic models
   │   │   ├── __init__.py
   │   │   ├── mean_field_solver.py
   │   │   ├── rootfinding_iterators.py
   │   │   ├── solver_base.py
   │   │   ├── steady_state_solver.py
   │   ├── tests                         # Unit test for SCAKS interfaces
   │   │   ├── __init__.py
   │   │   ├── correctors
   │   │   ├── mkm_inputs
   │   │   ├── models
   │   │   ├── parsers
   │   │   ├── plotters
   │   │   ├── solvers
   │   │   ├── test_all.py
   │   │   └── utilities
   │   └── utilities                     # Some utilities functions
   │       ├── __init__.py
   │       ├── check_utilities.py
   │       ├── coordinate_utilities.py
   │       ├── format_utilities.py
   │       └── profiling_utitlities.py
   ├── scripts                           # Some scripts for reference of users
   │   ├── mkm_plot                      # Some plotting scripts
   │   └── mkm_run                       # Run script
   ├── setup.cfg                         # Python setup configure file
   └── setup.py                          # Python setup script based on setuptools
