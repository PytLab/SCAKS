#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from scaks import __version__ as version

maintainer = 'Shao-Zheng-Jiang'
maintainer_email = 'shaozhengjiang@gmail.com'
author = maintainer
author_email = maintainer_email
description = "Micro-Kinetics Analysis package for Catalyst in Python"
long_description = '''
SCAKS
=====

Shanghai CAtlysis Kinetics Software

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

name = 'scaks'
platforms = ['linux', 'macos']
url = ''
download_url = ''

classifiers = [
    'Development Status :: 3 - Alpha',
    'Topic :: Utilities',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 2'
    'Programming Language :: Python :: 2.7'
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
]

test_suite = 'scaks.tests.test_all'

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
      version=version,
      test_suite=test_suite
  )

