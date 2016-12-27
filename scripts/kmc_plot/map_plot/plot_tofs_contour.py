#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

# Get pO2 dirs.
pO2_dirs = [i for i in commands.getoutput('ls').split('\n') if i.startswith('pO2-')]

pO2s = []
tofs = []
for pO2_dir in pO2_dirs:
    pO2 = float(pO2_dir.split('-')[-1])
    pO2s.append(pO2)
    # Get pCO dirs.
    cmd = "ls {}/".format(pO2_dir)
    pCO_dirs = [i for i in commands.getoutput(cmd).split('\n') if i.startswith('pCO-')]
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
        for idx in [0, 2, 8, 12]:
            TON += reaction_rates[reactions[idx]]
        tof = TON/3600
        tofs_1d.append(tof)
    tofs.append(tofs_1d)

pCOs = np.array(pCOs)
pO2s = np.array(pO2s)
tofs = np.array(tofs)

# 2D interpolate.
pO2s.shape = (-1, 1)
interp_func = interp2d(pCOs, pO2s, tofs, kind="linear")

# Plot 3D contour.
#y, x = np.mgrid[0:1:100j, 0:1:100j]
ynew = np.linspace(0.01, 2, 100)
xnew = np.linspace(0.01, 1, 100)
z = interp_func(xnew, ynew)

extent = [np.min(xnew), np.max(xnew), np.min(ynew), np.max(ynew)]

CS = plt.contourf(xnew.reshape(-1), ynew.reshape(-1),
                  z, 20, cmap=plt.cm.coolwarm)

CS2 = plt.contour(CS, levels=CS.levels[::2],
                  colors='#838B8B',
                  hold='on')
plt.clabel(CS2, colors='grey', inline=1, fontsize=12, fmt="%.2f")

# Add orthogonal lines.
# horizontal_line
#plt.plot([0.01, 1.0], [1.0, 1.0], color='#000000', linewidth=1.0)
#plt.plot([0.6, 0.6], [0.01, 2.0], color='#000000', linewidth=1.0)

plt.xlabel("P(CO_g)/bar")
plt.ylabel("P(O_2_g)/bar")

# Make a colorbar.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel("TOF/s^-1")

cbar.add_lines(CS2)

plt.show()

