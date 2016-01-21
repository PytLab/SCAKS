'''
    Module to plot ODE trajectory data in auto_ode_coverages.py
'''

import sys

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    adsorbate_names = sys.argv[1: ]
else:
    adsorbate_names = []

globs, locs = {}, {}
execfile('auto_ode_coverages.py', globs, locs)
times, coverages = locs['times'], locs['coverages']

# check args number
if adsorbate_names and len(adsorbate_names) != len(coverages[0]):
    print "args do not match coverges in auto_ode_coverages.py"
    sys.exit(1)

# plot
plt.subplot(111)
coverages = np.array(coverages)
if adsorbate_names:
    for idx in xrange(len(adsorbate_names)):
        ads = r'$\bf{' + adsorbate_names[idx] + r'^{*}}$'
        cvgs = coverages[:, idx]
        plt.plot(times, cvgs, label=ads, linewidth=2.5)
        plt.legend()
else:
    for idx in xrange(len(coverages[0])):
        cvgs = coverages[:, idx]
        plt.plot(times, cvgs, linewidth=2.5)

plt.ylabel(r'$\bf{Coverages}$', fontsize=16)
plt.xlabel(r'$\bf{Time (s)}$', fontsize=16)
plt.ylim(-0.1, 1.1)
plt.grid(True)

plt.show()

