#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

from kynetix.compatutil import subprocess

species_names = {"O": ("O_s", "O_r", "O_u", "O_d", "O_l"),
                 "CO": ("C")}

status, output = subprocess.getstatusoutput("find . -type d ! -name template")
dirs = output.split('\n')
dirs.remove('.')

pressures = []
coverages = {"O": [], "CO": []}

for d in dirs:
    filename = "{}/auto_coverages.py".format(d)
    if os.path.exists(filename):
        globs, locs = {}, {}
        exec(open(filename, "rb").read(), globs, locs)

        possible_types = locs["possible_types"]
        for species, type_names in species_names.items():
            coverage = 0.0
            for type_name in type_names:
                index = possible_types.index(type_name)
                coverage += locs["coverages"][index][-1]
            coverages[species].append(coverage)
        pressures.append(float(d[2:]))

p = np.array(pressures)
CO_cvgs = np.array(coverages["CO"])
O_cvgs = np.array(coverages["O"])

CO_cvgs = CO_cvgs[np.argsort(p)]
O_cvgs = O_cvgs[np.argsort(p)]
p = np.sort(p)

np.savez("cvgs-data.npz", p=p, CO_cvgs=CO_cvgs, O_cvgs=O_cvgs)

