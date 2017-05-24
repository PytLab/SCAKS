#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

with np.load("tofs-data.npz") as data:
    pCOs = data['pCOs']
    pO2s = data['pO2s']
    tofs = data['tofs']

with np.load("CO-cvgs-data.npz") as data:
    cvgs_CO = data['cvgs']

with np.load("O-cvgs-data.npz") as data:
    cvgs_O = data['cvgs']

cvgs = cvgs_O - cvgs_CO

# 2D interpolate.
pO2s.shape = (-1, 1)
interp_func_cvgs = interp2d(pCOs, pO2s, cvgs, kind="linear")
interp_func_tofs = interp2d(pCOs, pO2s, tofs, kind="linear")

# Plot 3D contour.
#y, x = np.mgrid[0:1:100j, 0:1:100j]
ynew = np.linspace(1e-5, 0.5, 100)
xnew = np.linspace(1e-5, 0.5, 200)
extent = [np.min(xnew), np.max(xnew), np.min(ynew), np.max(ynew)]
z_cvgs = interp_func_cvgs(xnew, ynew)
z_tofs = interp_func_tofs(xnew, ynew)

CS_cvgs = plt.contourf(xnew.reshape(-1), ynew.reshape(-1),
                       z_cvgs, 10, cmap=plt.cm.coolwarm)

CS_tofs = plt.contour(xnew.reshape(-1), ynew.reshape(-1), z_tofs, 10,
                      alpha=0.7, linewidths=1.0)

plt.xlabel(r"$P(CO_{(g)}/bar$")
plt.ylabel(r"$P(O_{2(g)}/bar$")

# Make a colorbar.
cbar = plt.colorbar(CS_tofs)
cbar.ax.set_ylabel("$TOF$")

#cbar.add_lines(CS2)

plt.show()

