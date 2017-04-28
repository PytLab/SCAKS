#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

from auto_cvgs import cvgs_O, cvgs_CO, pCO, pO2

pCO = np.array(pCO)
pO2 = np.array(pO2)
cvgs_O = np.array(cvgs_O)
cvgs_CO = np.array(cvgs_CO)

# 2D interpolate.
pO2.shape = (-1, 1)
interp_func = interp2d(pCO, pO2, cvgs_O, kind="linear")

# Plot 3D contour.
#y, x = np.mgrid[0:1:100j, 0:1:100j]
ynew = np.linspace(1e-5, 0.5, 100)
xnew = np.linspace(1e-5, 0.5, 100)
z = interp_func(xnew, ynew)

extent = [np.min(xnew), np.max(xnew), np.min(ynew), np.max(ynew)]

CS = plt.contourf(xnew.reshape(-1), ynew.reshape(-1),
                  z, 20, cmap=plt.cm.Reds)

CS2 = plt.contour(CS, levels=CS.levels[::2],
                  colors='#838B8B',
                  hold='on')
plt.clabel(CS2, colors='grey', inline=1, fontsize=12, fmt="%.2f")

plt.xlabel("P(CO_g)/bar")
plt.ylabel(r"P(O_2_g)/bar")

# Make a colorbar.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel("O_coverage")

cbar.add_lines(CS2)

plt.show()


