'''
    Script to plot auto_TOFs.py

    Usage example:
    --------------
    >>> python plot_TOFs.py auto_TOFs.py auto_TOFs_ctn.py
'''

import sys

import numpy as np
import matplotlib.pyplot as plt

# get filename
if len(sys.argv) < 2:
    filename = 'auto_TOFs.py' if len(sys.argv) == 1 else sys.argv[1]
    # get time and steps
    globs, locs = {}, {}
    execfile(filename, globs, locs)
    times, steps = locs['times'], locs['steps']

elif len(sys.argv) >= 2:
    globs, locs = {}, {}
    filename_1 = sys.argv[1]
    execfile(filename_1, globs, locs)

    # go through all the other files to append data
    for filename in sys.argv[2:]:
        temp_globs, temp_locs = {}, {}
        execfile(filename, temp_globs, temp_locs)
        # for different data
        for var in temp_locs:
            if var == 'times' or var == 'steps' or var.startswith('TOFs_'):
                locs[var].extend(temp_locs[var])
    times, steps = locs['times'], locs['steps']

else:
    print "Usage: python plot_TOFs.py file1 file2 ..."
    sys.exit(1)


# get TOFs
line_pts = []
log_line_pts = []
labels = []
log_labels = []
for var in locs:
    if var.startswith('TOFs_'):
        idx = var.split('_')[-1]
        forward_pts, reversed_pts = zip(*locs[var])
        log_forward_pts = np.log10(forward_pts)
        log_reversed_pts = np.log10(reversed_pts)

        # collect forward direction
        line_pts.append(forward_pts)
        log_line_pts.append(log_forward_pts)
        labels.append(idx + '^{+}')
        log_labels.append(idx + '^{+}')

        # collect reversed direction
        line_pts.append(reversed_pts)
        log_line_pts.append(log_reversed_pts)
        labels.append(idx + '^{-}')
        log_labels.append(idx + '^{-}')

if '--log' in sys.argv:
    # TOFs logrithm vs steps
    plt.subplot(121)
    for label, pts in zip(log_labels, log_line_pts):
        label = r'$\bf{' + label + r'}$'
        plt.plot(steps, pts, label=label, linewidth=2.5)

    plt.ylabel(r'$\bf{log(TOF)}$', fontsize=16)
    plt.xlabel(r'$\bf{kMC step}$', fontsize=16)
    plt.legend()
    plt.grid(True)

    # TOFs logrithm vs time
    plt.subplot(122)
    for label, pts in zip(log_labels, log_line_pts):
        label = r'$\bf{' + label + r'}$'
        plt.plot(times, pts, label=label, linewidth=2.5)

    plt.ylabel(r'$\bf{log(TOF)}$', fontsize=16)
    plt.xlabel(r'$\bf{Time (s)}$', fontsize=16)
    plt.legend()
    plt.grid(True)

else:
    # TOFs vs steps
    plt.subplot(221)
    for label, pts in zip(labels, line_pts):
        label = r'$\bf{' + label + r'}$'
        plt.plot(steps, pts, label=label, linewidth=2.5)

    plt.ylabel(r'$\bf{TOF}$', fontsize=16)
    plt.xlabel(r'$\bf{kMC step}$', fontsize=16)
    plt.legend()
    plt.grid(True)

    # TOFs vs time
    plt.subplot(222)
    for label, pts in zip(labels, line_pts):
        label = r'$\bf{' + label + r'}$'
        plt.plot(times, pts, label=label, linewidth=2.5)

    plt.ylabel(r'$\bf{TOF}$', fontsize=16)
    plt.xlabel(r'$\bf{Time (s)}$', fontsize=16)
    plt.legend()
    plt.grid(True)

    # TOFs logrithm vs steps
    plt.subplot(223)
    for label, pts in zip(log_labels, log_line_pts):
        label = r'$\bf{' + label + r'}$'
        plt.plot(steps, pts, label=label, linewidth=2.5)

    plt.ylabel(r'$\bf{log(TOF)}$', fontsize=16)
    plt.xlabel(r'$\bf{kMC step}$', fontsize=16)
    plt.legend()
    plt.grid(True)

    # TOFs logrithm vs time
    plt.subplot(224)
    for label, pts in zip(log_labels, log_line_pts):
        label = r'$\bf{' + label + r'}$'
        plt.plot(times, pts, label=label, linewidth=2.5)

    plt.ylabel(r'$\bf{log(TOF)}$', fontsize=16)
    plt.xlabel(r'$\bf{Time (s)}$', fontsize=16)
    plt.legend()
    plt.grid(True)

plt.show()
