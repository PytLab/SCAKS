#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt

with np.load("O-cvgs-data.npz") as data:
    pCOs = data['pCOs']
    pO2s = data['pO2s']
    cvgs = data['cvgs']

# 2D interpolate.
pO2s.shape = (-1, 1)
interp_func = interp2d(pCOs, pO2s, cvgs, kind="linear")

# Plot 3D contour.
#y, x = np.mgrid[0:1:100j, 0:1:100j]
ynew = np.linspace(1e-5, 0.5, 100)
xnew = np.linspace(1e-5, 0.5, 100)
z = interp_func(xnew, ynew)

extent = [np.min(xnew), np.max(xnew), np.min(ynew), np.max(ynew)]

CS = plt.contourf(xnew.reshape(-1), ynew.reshape(-1),
                  z, 10, cmap=plt.cm.Reds)

# Add orthogonal lines.
# horizontal_line
#plt.plot([0.01, 1.0], [1.0, 1.0], color='#000000', linewidth=1.0)
#plt.plot([0.6, 0.6], [0.01, 2.0], color='#000000', linewidth=1.0)

plt.xlabel("P(CO_g)/bar")
plt.ylabel("P(O_2_g)/bar")

# Make a colorbar.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel("O coverage")

plt.show()
#plt.savefig('tofs_contour.pdf')

