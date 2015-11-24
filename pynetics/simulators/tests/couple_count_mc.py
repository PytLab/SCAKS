import logging
import operator

import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import HexagonalSurface

ads = ('H', 'CO')
colors = ['#000000', '#F08080', '#228B22', '#4169E1',
          '#EE7621', '#BA55D3', '#708090', '#FFFF00', '#EE00EE']
co_cvg_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
line_list = []
h_cvgs_list = []
grid_length = 300

for co_cvg in co_cvg_list:
    h_cvgs = np.linspace(0, 1 - co_cvg, 71)
    h_cvgs_list.append(h_cvgs)
    probabilities = []
    for h_cvg in h_cvgs:
        surface = HexagonalSurface(grid_length, grid_length)
        surface.initialize_surface(ads, (h_cvg, co_cvg), ads)
        ncouple = surface.count_couples(ads)
        logging.info('ncouple = %d, h_cvg = %f', ncouple, h_cvg)
        #total_couple = operator.mul(*surface.shape)/2.0
        total_couple = grid_length**2
        p = float(ncouple)/total_couple
        logging.info('p = %f, h_cvg = %f', p, h_cvg)
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
