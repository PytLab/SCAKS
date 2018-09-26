#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from mikiac import __version__ as version

maintainer = 'Shao-Zheng-Jiang'
maintainer_email = 'shaozhengjiang@gmail.com'
author = maintainer
author_email = maintainer_email
description = "Micro-Kinetics Analysis package for Catalyst"
long_description = '''
MiKiAC
======

**Mi**\ cro-**Ki**\ netic **A**\ nalysis package for **C**\ atalyst

'''
install_requires = [
    'numpy',
    'scipy',
    'matplotlib',
    'sympy',
    'mpmath',
    'mpi4py',
    'prettytable'
]

license = 'LICENSE'

name = 'mikiac'
platforms = ['linux', 'windows', 'macos']
url = 'https://mikiac-docs.readthedocs.io/en/latest/'
download_url = 'https://github.com/PytLab/gaft/releases'

classifiers = [
    'Development Status :: 3 - Alpha',
    'Topic :: Utilities',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 2'
    'Programming Language :: Python :: 2.7'
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
]

setup(author=author,
      author_email=author_email,
      description=description,
      license=license,
      long_description=long_description,
      install_requires=install_requires,
      maintainer=maintainer,
      name=name,
      packages=find_packages(),
      platforms=platforms,
      url=url,
      download_url=download_url,
      version=version)

