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
exec(open('auto_ode_coverages.py', "rb").read(), globs, locs)
times, coverages = locs['times'], locs['coverages']

# check args number
if adsorbate_names and len(adsorbate_names) != len(coverages[0]):
    print("args do not match coverges in auto_ode_coverages.py")
    sys.exit(1)

# plot
fig = plt.figure()
ax = fig.add_subplot(111)

coverages = np.array(coverages)
CO_cvgs = coverages[:, 0]
O_cvgs = coverages[:, 1]
Vac_cvgs = np.ones(len(coverages)) - CO_cvgs - O_cvgs

ax.set_ylim(-0.1, 1.1)
time = times[-1]
left, right = -time*0.05, time*1.05
ax.set_xlim(left, right)

label = r'$\bf{CO^{*}}$'
ax.scatter(times, CO_cvgs, alpha=0.6, s=150.0, edgecolor="#000000", facecolor="#000000", label=label)
label = r'$\bf{O^{*}}$'
ax.scatter(times, O_cvgs, alpha=0.6, s=150.0, edgecolor="#E41615", facecolor="#E41615", label=label)
label = r'$\bf{*}$'
ax.scatter(times, Vac_cvgs, alpha=0.6, s=150.0, edgecolor="#0F387A", facecolor="#0F387A", label=label)
ax.plot([left, right], [1.916e-01, 1.916e-01], linewidth=1, color="#AEAEAE")
ax.plot([left, right], [3.322e-02, 3.322e-02], linewidth=1, color="#AEAEAE")

plt.legend()
plt.show()

