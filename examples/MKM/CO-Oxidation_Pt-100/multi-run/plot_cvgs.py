#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands
import os

import numpy as np
import matplotlib.pyplot as plt

colors = ["#CD6889", "#7D9EC0"]
species = {"CO": ["CO_s"], "O": ["O_s", "O2_2s"]}

if "__main__" == __name__:
    globs, locs = {}, {}
    execfile("auto_cvgs.py", globs, locs)
    pressures = locs['pCO']
    coverages = locs['cvgs']
    adsorbates = locs['adsorbates']


    fig = plt.figure()
    ax = fig.add_subplot(111)
    pressure_range = max(pressures) - min(pressures)
    ax.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel(r"$\bf{P(CO_g)/bar}$")
    ax.set_ylabel(r"$\bf{Coverages}$")


    coverages = zip(*coverages)
    CO_cvgs = np.array(coverages[0]) + np.array(coverages[1])
    O_cvgs = np.array(coverages[2])
    cvgs_dict = {"CO": CO_cvgs, "O": O_cvgs}
    for idx, sp in enumerate(cvgs_dict):
        ax.plot(pressures, cvgs_dict[sp], color=colors[idx],
                                     linewidth=2.0,
                                     marker='o',
                                     markerfacecolor=colors[idx],
                                     markeredgecolor=colors[idx],
                                     markersize=7.0,
                                     label=sp)
    ax.legend()
    ax.grid(True)
    plt.show()

