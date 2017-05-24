#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

with np.load("tofs-data.npz") as data:
    pressures = data['p']
    tofs = data['tofs']

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
pressure_range = max(pressures) - min(pressures)
ax.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
tof_range = max(tofs) - min(tofs)
ax.set_ylim(min(tofs)-tof_range*0.1, max(tofs)+tof_range*0.1)
ax.set_xlabel(r"$\bf{P(CO_g)/bar}$")
ax.set_ylabel(r"$\bf{TOF/s^-1}$")
ax.plot(pressures, tofs, color="#7D9EC0",
                         linewidth=2.0,
                         marker='o',
                         markerfacecolor="#CD6889",
                         markeredgecolor="#CD6889",
                         markersize=7.0)
ax.grid(True)
plt.show()

