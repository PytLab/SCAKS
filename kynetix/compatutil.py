#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for wrap functions for compatibility of Python2 and Python3
"""
import sys

# Check version for current python interpreter.
if sys.version > "3":
    PY2 = False
else:
    PY2 = True

# Compatible reduce function.
if PY2:
    reduce = reduce
else:
    from functools import reduce

# Compatible Exception class.
if PY2:
    from exceptions import Exception
else:
    Exception = Exception

