import matplotlib.pyplot as plt

globs, locs = {}, {}
execfile('auto_coverages.py', globs, locs)

times, steps, coverages = locs['times'], locs['steps'], locs['coverages']
possible_types = locs['possible_types']

# coverage vs steps
plt.subplot(121)
for sp, cvgs in zip(possible_types, coverages):
    plt.plot(steps, cvgs, label=sp, linewidth=2.5)

plt.ylabel(r'$Coverages$')
plt.xlabel(r'$kMC step$')
plt.ylim(-0.1, 1.1)

plt.subplot(122)
for sp, cvgs in zip(possible_types, coverages):
    plt.plot(times, cvgs, label=sp, linewidth=2.5)

plt.ylabel(r'$Coverages$')
plt.xlabel(r'$Time (s)$')

plt.ylim(-0.1, 1.1)
plt.legend()
plt.show()
