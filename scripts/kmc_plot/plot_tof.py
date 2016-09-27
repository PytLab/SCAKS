'''
Module to plot auto_tofs.py
'''

import sys

import matplotlib.pyplot as plt

from kynetix.parsers.rxn_parser import RxnEquation

production = "CO2_g"

if "__main__" == __name__:
    globs, locs = {}, {}
    execfile("auto_tofs.py", globs, locs)
    rates_lists = locs["tofs"]
    times = locs["times"]
    processes = locs["processes"]

    # Get forward and reverse tof indices.
    forward_indices, reverse_indices = [], []
    for idx, process in enumerate(processes):
        # Forward.
        if process.endswith("(->)"):
            rxn_str = process.strip("(->)")
            rxn_equation = RxnEquation(rxn_str)
            final_state = rxn_equation.tolist()[-1].chem_state()
            if production in final_state:
                forward_indices.append(idx)
        elif process.endswith("(<-)"):
            rxn_str = process.strip("(<-)")
            rxn_equation = RxnEquation(rxn_str)
            final_state = rxn_equation.tolist()[-1].chem_state()
            if production in final_state:
                reverse_indices.append(idx)

    # Get TOF.
    tofs = []
    for rates in rates_lists:
        tof = 0.0 
        for idx in forward_indices:
            tof += rates[idx]
        for idx in reverse_indices:
            tof -= rates[idx]
        tofs.append(tof)

    plt.subplot(111)
    plt.plot(times, tofs, linewidth=2.5)

    plt.ylabel(r"$\bf{TOF}$", fontsize=16)
    plt.xlabel(r"$\bf{time}$", fontsize=16)
    plt.legend()
    plt.grid(True)

    if "-s" in sys.argv:
        plt.savefig("tofs.png")
    else:
        plt.show()

