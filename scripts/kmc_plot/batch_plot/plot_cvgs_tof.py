#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

with np.load("tofs-data.npz") as data:
    pressures = data['p']
    tofs = data['tofs']

with np.load("cvgs-data.npz") as data:
    CO_cvgs = data['CO_cvgs']
    O_cvgs = data['O_cvgs']

fig = plt.figure()
ax = fig.add_subplot(111)
pressure_range = max(pressures) - min(pressures)
ax.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
ax.set_ylim(-0.1, 1.1)
ax.set_xlabel(r"$P_{CO_{(g)}}/bar$")
ax.set_ylabel(r"$Coverages$")

# Plot CO
ax.plot(pressures, CO_cvgs, linewidth=2.0,
                            color="#0F0F0F",
                            marker='o',
                            markerfacecolor="#0F0F0F",
                            markeredgecolor="#0F0F0F",
                            markersize=7.0,
                            alpha=0.7,
                            label=r"$\theta_{CO^{*}}$")

# Plot O
ax.plot(pressures, O_cvgs, linewidth=2.0,
                           color="#CD6889",
                           marker='o',
                           markerfacecolor="#CD6889",
                           markeredgecolor="#CD6889",
                           markersize=7.0,
                           alpha=0.7,
                           label=r"$\theta_{O^{*}}$")

# Plot TOF.
ax2 = ax.twinx()
tof_range = max(tofs) - min(tofs)
ax2.set_ylim(min(tofs)-tof_range*0.1, max(tofs)+tof_range*0.1)
ax2.set_ylabel(r"$TOF/s^{-1}$")
ax2.plot(pressures, tofs, color="#7D9EC0",
                          linewidth=2.0,
                          marker='s',
                          markerfacecolor="#7D9EC0",
                          markeredgecolor="#7D9EC0",
                          markersize=7.0,
                          alpha=0.3,
                          label="$TOF$")
# Trick for labels.
ax.plot(0, 0, color="#0F0F0F",
              linewidth=2.0,
              marker='s',
              markerfacecolor="#0F0F0F",
              markeredgecolor="#0F0F0F",
              markersize=7.0,
              alpha=0.3,
              label='$TOF$')
ax.legend(loc=0)  
#ax.grid(True)

plt.show()

