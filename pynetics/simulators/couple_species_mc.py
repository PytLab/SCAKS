import logging

import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import HexagonalSurface, get_probabilities

#logging setting
logging.basicConfig(format='%(name)s  %(levelname)s  %(message)s', level=logging.INFO)

ads = ('H', 'CO')
co_cvg = 0.5
h_cvgs = np.linspace(0, 1 - co_cvg, 71)
colors = ['#000000', '#F08080', '#228B22', '#4169E1',
          '#EE7621', '#BA55D3', '#708090']
target_numbers = range(7)
probabilities_list = [[]]*7

for h_cvg in h_cvgs:
    surface = HexagonalSurface(100, 100)
    surface.initialize_surface(ads, (h_cvg, co_cvg), ads)
    probabilities = get_probabilities(surface, 'H', 'H', target_numbers)
    logging.info('coverage = %.3f, probabilities = %s', h_cvg, probabilities)
    # append colums(axis=1)
    probabilities = np.array(probabilities).reshape(-1, 1)
    probabilities_list = np.append(probabilities_list, probabilities, axis=1)

# plot
for idx, probabilities in enumerate(probabilities_list):
    plt.scatter(h_cvgs, probabilities, s=20, color=colors[idx],
                alpha=0.8, label=('nn = %d' % idx))
s_co_cvg = '%.2f' % co_cvg
plt.xlabel(r'$\theta_{H}  (with \theta_{CO}=' + s_co_cvg + ')$')
plt.ylabel(r'$Probability$')
plt.grid(True)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.2)
plt.legend()

plt.show()
