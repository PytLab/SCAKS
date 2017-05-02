import os

import numpy as np
import matplotlib.pyplot as plt


if "__main__" == __name__:
    # Plot TOFs.
    globs, locs = {}, {}
    exec(open("auto_tofs.py", "rb").read(), globs, locs)
    pressures = locs['pO2']
    TONs = locs['tofs']

    fig = plt.figure()

    # Plot cvgs.
    globs, locs = {}, {}
    exec(open("auto_cvgs.py", "rb").read(), globs, locs)
    coverages = locs['cvgs']
    adsorbates = locs['adsorbates']


    ax1 = fig.add_subplot(111)
    pressure_range = max(pressures) - min(pressures)
    ax1.set_ylim(-0.1, 1.1)
    ax1.set_xlim(min(pressures)-pressure_range*0.1, max(pressures)+pressure_range*0.1)
    ax1.set_xlabel(r"$\bf{P(CO_g)/bar}$")
    ax1.set_ylabel(r"$\bf{Coverages}$")


    coverages = list(zip(*coverages))
    CO_cvgs = np.array(coverages[0])
    O_cvgs = np.array(coverages[1])
    cvgs_dict = {"CO": CO_cvgs, "O": O_cvgs}
    colors = ["#CD6889", "#7D9EC0"]
    for idx, sp in enumerate(cvgs_dict):
        ax1.plot(pressures, cvgs_dict[sp], color=colors[idx],
                                           linewidth=2.0,
                                           marker='o',
                                           markerfacecolor=colors[idx],
                                           markeredgecolor=colors[idx],
                                           markersize=14.0,
                                           alpha=0.7,
                                           label=sp)

    # Plot TOF.
    ax2 = ax1.twinx()
    pressure_range = max(pressures) - min(pressures)
    TON_range = max(TONs) - min(TONs)
    ax2.set_ylim(min(TONs)-TON_range*0.1, max(TONs)+TON_range*0.1)
    ax2.set_ylabel(r"$\bf{TOF/s^-1}$")
    ax2.plot(pressures, TONs, color="#0F0F0F",
                              linewidth=2.0,
                              marker='^',
                              markerfacecolor="#0F0F0F",
                              markeredgecolor="#0F0F0F",
                              markersize=14.0,
                              alpha=0.3,
                              label='TOF')

    # Trick for labels.
    ax1.plot(0, 0, color="#0F0F0F",
                   linewidth=2.0,
                   marker='^',
                   markerfacecolor="#0F0F0F",
                   markeredgecolor="#0F0F0F",
                   markersize=14.0,
                   alpha=0.3,
                   label='TOF')
    ax1.grid(True)
    ax1.legend(loc=0)

    plt.show()

