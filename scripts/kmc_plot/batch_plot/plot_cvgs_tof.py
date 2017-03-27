#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands
import os

import matplotlib.pyplot as plt

if "__main__" == __name__:
    status, output = commands.getstatusoutput("ls")
    dirs = [d for d in output.split("\n") if d.startswith("0.")]

    # Plot cvgs.
    species_names = {"O": ("O_s", "O_r", "O_u", "O_d", "O_l"),
                     "CO": ("C")}
    colors = ["#CD6889", "#7D9EC0"]

    pressures = []
    coverages = {"O": [], "CO": []}

    for d in dirs:
        filename = "{}/auto_coverages.py".format(d)
        if os.path.exists(filename):
            globs, locs = {}, {}
            exec(open(filename, "rb").read(), globs, locs)

            possible_types = locs["possible_types"]
            for species, type_names in species_names.items():
                coverage = 0.0
                for type_name in type_names:
                    index = possible_types.index(type_name)
                    coverage += locs["coverages"][index][-1]
                coverages[species].append(coverage)
            pressures.append(float(d))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    pressure_range = max(pressures) - min(pressures)
    ax.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel(r"$\bf{P(CO_g)/bar}$")
    ax.set_ylabel(r"$\bf{Coverages}$")
    for idx, species in enumerate(coverages):
        ax.plot(pressures, coverages[species], linewidth=2.0,
                                               color=colors[idx],
                                               marker='o',
                                               markerfacecolor=colors[idx],
                                               markeredgecolor=colors[idx],
                                               markersize=14.0,
                                               label=species,
                                               alpha=0.7)

    # Plot TOF.
    TONs = []
    for d in dirs:
        filename = "{}/auto_frequency.py".format(d)
        if os.path.exists(filename):
            globs, locs = {}, {}
            exec(open(filename, "rb").read(), globs, locs)
            try:
                reaction_rates = locs["reaction_rates"]
            except KeyError:
                continue
            TON = 0.0
            reactions = sorted(reaction_rates.keys())
            for idx in [0, 2, 8, 12]:
                TON += reaction_rates[reactions[idx]]
            TONs.append(TON)

    # Plot
    ax2 = ax.twinx()
    TON_range = max(TONs) - min(TONs)
    ax2.set_ylim(min(TONs)-TON_range*0.1, max(TONs)+TON_range*0.1)
    ax2.set_ylabel(r"$\bf{TON/s^-1}$")
    ax2.plot(pressures, TONs, color="#0F0F0F",
                              linewidth=2.0,
                              marker='^',
                              markerfacecolor="#0F0F0F",
                              markeredgecolor="#0F0F0F",
                              markersize=14.0,
                              alpha=0.3,
                              label="TOF")
    # Trick for labels.
    ax.plot(0, 0, color="#0F0F0F",
                  linewidth=2.0,
                  marker='^',
                  markerfacecolor="#0F0F0F",
                  markeredgecolor="#0F0F0F",
                  markersize=14.0,
                  alpha=0.3,
                  label='TOF')
    ax.legend(loc=0)  
    ax.grid(True)

    plt.show()

