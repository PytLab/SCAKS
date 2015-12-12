import logging
import operator

import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import SquareSurface

ads = ('H', )
colors = ['#000000', '#F08080', '#228B22', '#4169E1',
          '#EE7621', '#BA55D3', '#708090', '#FFFF00', '#EE00EE']
h_cvgs = np.linspace(0.0, 1.0, 100)

probabilities = []
for h_cvg in h_cvgs:
#    surface = SquareSurface(200, 200)
#    surface.initialize_surface(ads, (h_cvg, ), ads)
#    ncouple = surface.count_couples(ads)
#    logging.info('ncouple = %d, h_cvg = %f', ncouple, h_cvg)
#    total_couple = operator.mul(*surface.shape)/2.0
#    p = float(ncouple)/total_couple
    p = h_cvg*(1 - (1 - h_cvg)**2)
    logging.info('p = %f, h_cvg = %f', p, h_cvg)
    probabilities.append(p)

plt.scatter(h_cvgs, probabilities, s=10, color='#BA55D3', alpha=0.8)

plt.xlabel(r'$\bf{\theta_{H}}$')
plt.ylabel('Probability of finding couples')
plt.grid(True)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.2)
plt.legend()

plt.show()
