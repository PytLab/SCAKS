'''
    Module to plot auto_converages.py
'''
import sys

import numpy as np
import matplotlib.pyplot as plt

globs, locs = {}, {}
exec(open('auto_frequency.py', "rb").read(), globs, locs)

reaction_occurencies = locs['reaction_occurencies']
reactions = sorted(reaction_occurencies)

# Get data.
N = len(reaction_occurencies)/2
forward_rxns = reactions[::2]
reverse_rxns = reactions[1::2]
f_occurencies = [reaction_occurencies[rxn] for rxn in forward_rxns]
r_occurencies = [reaction_occurencies[rxn] for rxn in reverse_rxns]

ind = np.arange(N)
width = 0.5

fig = plt.figure()
ax = fig.add_subplot(111)

rects1 = ax.barh(ind, f_occurencies, width, color='#CD8C95', linewidth=0)
rects2 = ax.barh(ind + width, r_occurencies, width, color='#7D9EC0', linewidth=0)

ax.set_ylim(-0.5, N + width + 0.5)
plt.xlabel(r'$\bf{Frequency}$', fontsize=16)
ax.set_yticks(ind + width)

# Get yticklabels.
ylabels = [r"$\bf{reaction(%d)}$" % i for i in xrange(N)]
ax.set_yticklabels(ylabels, rotation=-10, fontsize=13)

plt.legend((rects1[0], rects2[0]), (r'$\bf{forward}$', r'$\bf{reverse}$'))

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        width = rect.get_width()
        height = rect.get_height()
        ax.text(width + 0.3, rect.get_y() + height/4., r'$\bf{%d}$' % int(width), fontsize=13)

autolabel(rects1)
autolabel(rects2)

plt.grid(True)

if "-s" in sys.argv:
    plt.savefig("frequencies.png")
else:
    plt.show()

