#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands
import os

import numpy as np
import matplotlib.pyplot as plt

from scaks.compatutil import subprocess

status, output = subprocess.getstatusoutput("find . -type d ! -name template")
dirs = output.split('\n')
dirs.remove('.')

tofs = []
pressures = []
for d in dirs:
    filename = "{}/auto_frequency.py".format(d)
    if os.path.exists(filename):
        globs, locs = {}, {}
        exec(open(filename, "rb").read(), globs, locs)
        try:
            reaction_rates = locs["reaction_rates"]
        except KeyError:
            continue
        tof = 0.0
        reactions = sorted(reaction_rates.keys())
        for idx in [0, 2]:
            tof += reaction_rates[reactions[idx]]
        tofs.append(tof)
        pressures.append(float(d[2:]))

p = np.array(pressures)
tofs = np.array(tofs)

# Sort.
tofs = tofs[np.argsort(p)]
p = np.sort(p)

np.savez("tofs-data.npz", p=p, tofs=tofs)

