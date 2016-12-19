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
interp_func = interp2d(pCO, pO2, tofs, kind="cubic")

# Plot 3D contour.
y, x = np.ogrid[0:1:100j, 0:1:100j]
z = interp_func(x, y)

extent = [np.min(x), np.max(x), np.min(y), np.max(y)]

plt.contourf(x.reshape(-1), y.reshape(-1), z, 20)
plot.show()

