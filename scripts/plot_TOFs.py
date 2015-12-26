'''
    Module to plot auto_TOFs.py
'''

import matplotlib.pyplot as plt

globs, locs = {}, {}
execfile('auto_TOFs.py', globs, locs)

# time and steps
times, steps = locs['times'], locs['steps']
# get TOFs
line_pts = []
labels = []
for var in locs:
    if var.startswith('TOFs_'):
        idx = var.split('_')[-1]
        forward_pts, reversed_pts = zip(*locs[var])
        # collect forward direction
        line_pts.append(forward_pts)
        labels.append(idx + '^{+}')
        # collect reversed direction
        line_pts.append(reversed_pts)
        labels.append(idx + '^{-}')

# TOFs vs steps
plt.subplot(121)
for label, pts in zip(labels, line_pts):
    label = r'$\bf{' + label + r'}$'
    plt.plot(steps, pts, label=label, linewidth=2.5)

plt.ylabel(r'$\bf{TOFs}$', fontsize=16)
plt.xlabel(r'$\bf{kMC step}$', fontsize=16)
plt.legend()
plt.grid(True)

# TOFs vs time
plt.subplot(122)
for label, pts in zip(labels, line_pts):
    label = r'$\bf{' + label + r'}$'
    plt.plot(times, pts, label=label, linewidth=2.5)

plt.ylabel(r'$\bf{TOFs}$', fontsize=16)
plt.xlabel(r'$\bf{Time (s)}$', fontsize=16)
plt.legend()
plt.grid(True)

plt.show()
