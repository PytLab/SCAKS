'''
bond concentration MC simulation.
'''
import logging

import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import HexagonalSurface, count_neighbors

ads = ('H', 'CO')
colors = ['#000000', '#F08080', '#228B22', '#4169E1',
          '#EE7621', '#BA55D3', '#708090', '#FFFF00', '#EE00EE']
co_cvg_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
line_list = []
h_cvgs_list = []
tot_bond = 200*200*6/2.0

for co_cvg in co_cvg_list:
    h_cvgs = np.linspace(0, 1 - co_cvg, 71)
    h_cvgs_list.append(h_cvgs)
    probabilities = []
    for h_cvg in h_cvgs:
        logging.info('H cvg = %.2f, CO cvg = %.2f', h_cvg, co_cvg)
        surface = HexagonalSurface(200, 200)
        surface.initialize_surface(ads, (h_cvg, co_cvg), ads)
        nbond = 0
        for row in xrange(200):
            for col in xrange(200):
                site = surface.grid[row][col]
                if site and site.name == 'H':
                    nbond += count_neighbors(surface, row, col, 'H')
        nbond /= 2.0
        p = nbond/tot_bond
        logging.info('H bond concentration = %.2f, H concentration = %.2f', p, h_cvg)
        probabilities.append(p)
    line_list.append(probabilities)

for i in xrange(9):
    plt.scatter(h_cvgs_list[i], line_list[i], s=10, color=colors[i],
                alpha=0.8, label=(r'$\theta_{CO} = %.2f$' % co_cvg_list[i]))

plt.xlabel(r'$\bf{\theta_{H}}$')
plt.ylabel('Probability of finding couples')
plt.grid(True)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.2)
plt.legend()

plt.show()
