#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands
import matplotlib.pyplot as plt

if "__main__" == __name__:
    status, output = commands.getstatusoutput("ls")
    dirs = [d for d in output.split("\n") if d.startswith("1.")]

    TONs = []
    pressures = []
    for d in dirs:
        pressures.append(float(d) - 1)
        filename = "{}/auto_frequency.py".format(d)
        globs, locs = {}, {}
        execfile(filename, globs, locs)
        reaction_rates = locs["reaction_rates"]
        TON = 0.0
        reactions = sorted(reaction_rates.keys())
        for idx in [0, 2, 8, 12]:
            TON += reaction_rates[reactions[idx]]
        TONs.append(TON)

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pressure_range = max(pressures) - min(pressures)
    ax.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
    TON_range = max(TONs) - min(TONs)
    ax.set_ylim(min(TONs)-TON_range*0.1, max(TONs)+TON_range*0.1)
    ax.set_xlabel(r"$\bf{P(CO_g)/bar}$")
    ax.set_ylabel(r"$\bf{TON/s^-1}$")
    ax.plot(pressures, TONs, color="#7D9EC0",
                             linewidth=2.0,
                             marker='o',
                             markerfacecolor="#CD6889",
                             markeredgecolor="#CD6889",
                             markersize=7.0)
    ax.grid(True)
    plt.show()

