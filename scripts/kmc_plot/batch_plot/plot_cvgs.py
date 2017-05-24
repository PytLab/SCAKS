#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

with np.load("cvgs-data.npz") as data:
    pressures = data["p"]
    CO_cvgs = data['CO_cvgs']
    O_cvgs = data['O_cvgs']

# Plot.
fig = plt.figure()
ax = fig.add_subplot(111)
pressure_range = max(pressures) - min(pressures)
ax.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
ax.set_ylim(-0.1, 1.1)
ax.set_xlabel(r"$P_{CO_{(g)}}/bar$")
ax.set_ylabel(r"$Coverages$")

# Plot CO cvgs.
ax.plot(pressures, CO_cvgs, linewidth=2.0,
                            color="#CD6889",
                            marker='o',
                            markerfacecolor="#CD6889",
                            markeredgecolor="#CD6889",
                            markersize=7.0,
                            label=r"$\theta_{CO^{*}}$")

# Plot O cvgs.
ax.plot(pressures, O_cvgs, linewidth=2.0,
                           color="#7D9EC0",
                           marker='o',
                           markerfacecolor="#7D9EC0",
                           markeredgecolor="#7D9EC0",
                           markersize=7.0,
                           label=r"$\theta_{O^{*}}$")
ax.legend()  
plt.show()

