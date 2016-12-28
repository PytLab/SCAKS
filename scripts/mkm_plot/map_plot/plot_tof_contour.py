#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

from auto_tofs import pCO, pO2, tofs

pCO = np.array(pCO)
pO2 = np.array(pO2)
tofs = np.array(tofs)

# 2D interpolate.
pO2.shape = (-1, 1)
interp_func = interp2d(pCO, pO2, tofs, kind="linear")

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
plt.plot([0.01, 1.0], [1.0, 1.0], color='#000000', linewidth=1.0)
plt.plot([0.6, 0.6], [0.01, 2.0], color='#000000', linewidth=1.0)

plt.xlabel("P(CO_g)/bar")
plt.ylabel("P(O_2_g)/bar")

# Make a colorbar.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel("TOF/s^-1")

cbar.add_lines(CS2)

plt.show()

