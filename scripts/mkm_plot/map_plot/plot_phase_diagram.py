#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

from auto_cvgs import cvgs_O, cvgs_CO, pCO, pO2
from auto_tofs import tofs

pCO = np.array(pCO)
pO2 = np.array(pO2)
cvgs_O = np.array(cvgs_O)
cvgs_CO = np.array(cvgs_CO)
cvgs = cvgs_O - cvgs_CO
tofs = np.array(tofs)

# 2D interpolate.
pO2.shape = (-1, 1)
interp_func_cvgs = interp2d(pCO, pO2, cvgs, kind="linear")
interp_func_tofs = interp2d(pCO, pO2, tofs, kind="linear")

# Plot 3D contour.
#y, x = np.mgrid[0:1:100j, 0:1:100j]
ynew = np.linspace(1e-5, 0.5, 100)
xnew = np.linspace(1e-5, 0.5, 100)
extent = [np.min(xnew), np.max(xnew), np.min(ynew), np.max(ynew)]
z_cvgs = interp_func_cvgs(xnew, ynew)
z_tofs = interp_func_tofs(xnew, ynew)

CS_cvgs = plt.contourf(xnew.reshape(-1), ynew.reshape(-1),
                      z_cvgs, 5, cmap=plt.cm.coolwarm)

CS_tofs = plt.contour(xnew.reshape(-1), ynew.reshape(-1), z_tofs, 10,
                      alpha=0.7, linewidths=1.0)

plt.xlabel(r"$P(CO_{(g)})/bar$")
plt.ylabel(r"$P(O_{2(g)}/bar$")

# Make a colorbar.
cbar = plt.colorbar(CS_tofs)
cbar.ax.set_ylabel(r"$TOF$")

#cbar.add_lines(CS2)
# Add lines.
plt.plot([0.1, 0.1], [0.2, 0.5], lw=1, color="#FFFFFF", marker="_")
plt.plot([0.0, 0.2], [0.35, 0.35], lw=1, color="#C5C5C2", marker="|")

# set limits.
plt.xlim([0.0, 0.5])
plt.ylim([0.0, 0.5])


plt.show()


