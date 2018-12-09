#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from scaks.compatutil import subprocess

# Get pO2 dirs.
pO2_dirs = [i for i in subprocess.getoutput('ls').split('\n') if i.startswith('pO2-')]

pO2s = []
cvgs = []
for pO2_dir in pO2_dirs:
    pO2 = float(pO2_dir.split('-')[-1])
    pO2s.append(pO2)
    # Get pCO dirs.
    cmd = "ls {}/".format(pO2_dir)
    pCO_dirs = [i for i in subprocess.getoutput(cmd).split('\n') if i.startswith('pCO-')]
    pCOs = [float(pCO_dir.split('-')[-1]) for pCO_dir in pCO_dirs]
    cvgs_1d = []
    for pCO_dir in pCO_dirs:
        # Read tofs.
        filename = "{}/{}/auto_coverages.py".format(pO2_dir, pCO_dir)
        globs, locs = {}, {}
        execfile(filename, globs, locs)
        cvg = locs["coverages"][-1][-1]
        cvgs_1d.append(cvg)
    cvgs.append(cvgs_1d)

pCOs = np.array(pCOs)
pO2s = np.array(pO2s)
cvgs = np.array(cvgs)

np.savez("CO-cvgs-data.npz", pCOs=pCOs, pO2s=pO2s, cvgs=cvgs)

