#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import matplotlib.pyplot as plt

if "__main__" == __name__:
    globs, locs = {}, {}
    exec(open("auto_tofs.py", "rb").read(), globs, locs)
    pressures = locs['pCO']
    TONs = locs['tofs']

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pressure_range = max(pressures) - min(pressures)
    ax.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
    TON_range = max(TONs) - min(TONs)
    ax.set_ylim(min(TONs)-TON_range*0.1, max(TONs)+TON_range*0.1)
    ax.set_xlabel(r"$\bf{P(CO_g)/bar}$")
    ax.set_ylabel(r"$\bf{TOF/s^-1}$")
    ax.plot(pressures, TONs, color="#7D9EC0",
                             linewidth=2.0,
                             marker='o',
                             markerfacecolor="#CD6889",
                             markeredgecolor="#CD6889",
                             markersize=7.0)
    ax.grid(True)
    plt.show()

