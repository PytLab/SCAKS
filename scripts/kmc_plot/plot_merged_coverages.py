'''
    Module to plot auto_converages.py
'''

import sys

import numpy as np
import matplotlib.pyplot as plt

if "__main__" == __name__:
    globs, locs = {}, {}
    exec(open('auto_coverages.py', "rb").read(), globs, locs)

    times, steps, coverages = locs['times'], locs['steps'], locs['coverages']
    possible_types = locs['possible_types']

    # coverage vs steps
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    CO_cvgs, O_cvgs, Vac_cvgs = np.zeros(len(steps)), np.zeros(len(steps)), np.zeros(len(steps))
    for sp, cvgs in zip(possible_types, coverages):
        if sp.startswith("O"):
            O_cvgs += np.array(cvgs)
        elif sp.startswith("C"):
            CO_cvgs += np.array(cvgs)
        elif sp.startswith("V"):
            Vac_cvgs += np.array(cvgs)

    O_avg = np.mean(O_cvgs)
    CO_avg = np.mean(CO_cvgs)

    # Plot lines.
    # CO
    label = r'$\bf{CO^{*}}$'
    ax1.scatter(steps, CO_cvgs, alpha=0.3, facecolor="#93989A", edgecolor="#000000", label=label)
    # O
    label = r'$\bf{O^{*}}$'
    ax1.scatter(steps, O_cvgs, alpha=0.3, facecolor="#E0E0E0", edgecolor="#E41615", label=label)
    # Vac
    label = r'$\bf{*}$'
    ax1.scatter(steps, Vac_cvgs, alpha=0.3, facecolor="#93989A", edgecolor="#0F387A", label=label)

    ax1.set_ylim(-0.1, 1.1)
    step = steps[-1]
    left, right = -step*0.05, step*1.05
    ax1.set_xlim(left, right)

    # average lines.
    ax1.plot([left, right], [CO_avg, CO_avg], linewidth=1, color="#AEAEAE")
    ax1.plot([left, right], [O_avg, O_avg], linewidth=1, color="#AEAEAE")

    # coverage vs time
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    # Plot lines.
    # CO
    label = r'$\bf{CO^{*}}$'
    ax2.scatter(times, CO_cvgs, alpha=0.3, facecolor="#93989A", edgecolor="#000000", label=label)
    # O
    label = r'$\bf{O^{*}}$'
    ax2.scatter(times, O_cvgs, alpha=0.3, facecolor="#E0E0E0", edgecolor="#E41615", label=label)
    # Vac
    label = r'$\bf{*}$'
    ax2.scatter(times, Vac_cvgs, alpha=0.3, facecolor="#93989A", edgecolor="#0F387A", label=label)

    ax2.set_ylim(-0.1, 1.1)
    time = times[-1]
    left, right = -time*0.05, time*1.05
    ax2.set_xlim(left, right)
    ax2.plot([left, right], [CO_avg, CO_avg], linewidth=1, color="#AEAEAE")
    ax2.plot([left, right], [O_avg, O_avg], linewidth=1, color="#AEAEAE")

    plt.show()

