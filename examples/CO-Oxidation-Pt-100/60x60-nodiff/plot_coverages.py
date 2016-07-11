'''
    Module to plot auto_converages.py
'''

import matplotlib.pyplot as plt

globs, locs = {}, {}
execfile('auto_coverages.py', globs, locs)

times, steps, coverages = locs['times'], locs['steps'], locs['coverages']
possible_types = locs['possible_types']

# coverage vs steps
plt.subplot(121)
for sp, cvgs in zip(possible_types, coverages):
    sp = r'$\bf{*}$' if sp == 'Vac' else r'$\bf{' + sp + r'^{*}}$'
    plt.plot(steps, cvgs, label=sp, linewidth=2.5)

plt.ylabel(r'$\bf{Coverages}$', fontsize=16)
plt.xlabel(r'$\bf{kMC step}$', fontsize=16)
plt.ylim(-0.1, 1.1)
plt.legend()
plt.grid(True)

# coverage vs time
plt.subplot(122)
for sp, cvgs in zip(possible_types, coverages):
    sp = r'$\bf{*}$' if sp == 'Vac' else r'$\bf{' + sp + r'^{*}}$'
    plt.plot(times, cvgs, label=sp, linewidth=2.5)

plt.ylabel(r'$\bf{Coverages}$', fontsize=16)
plt.xlabel(r'$\bf{Time (s)}$', fontsize=16)
plt.ylim(-0.1, 1.1)
plt.legend()
plt.grid(True)

plt.show()
