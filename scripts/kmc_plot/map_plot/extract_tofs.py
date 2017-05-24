#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

from kynetix.compatutil import subprocess

# Get pO2 dirs.
pO2_dirs = [i for i in subprocess.getoutput('ls').split('\n') if i.startswith('pO2-')]

pO2s = []
tofs = []
for pO2_dir in pO2_dirs:
    pO2 = float(pO2_dir.split('-')[-1])
    pO2s.append(pO2)
    # Get pCO dirs.
    cmd = "ls {}/".format(pO2_dir)
    pCO_dirs = [i for i in subprocess.getoutput(cmd).split('\n') if i.startswith('pCO-')]
    pCOs = [float(pCO_dir.split('-')[-1]) for pCO_dir in pCO_dirs]
    tofs_1d = []
    for pCO_dir in pCO_dirs:
        # Read tofs.
        filename = "{}/{}/auto_frequency.py".format(pO2_dir, pCO_dir)
        globs, locs = {}, {}
        execfile(filename, globs, locs)
        reaction_rates = locs["reaction_rates"]
        TON = 0.0
        reactions = sorted(reaction_rates.keys())
        for idx in [0, 2, 8]:
            TON += reaction_rates[reactions[idx]]
        tof = TON
        tofs_1d.append(tof)
    tofs.append(tofs_1d)

pCOs = np.array(pCOs)
pO2s = np.array(pO2s)
tofs = np.array(tofs)

# Write data to file.
np.savez("tof-data.npz", pCOs=pCOs, pO2s=pO2s, tofs=tofs)

