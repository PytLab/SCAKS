import logging

import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import HexagonalSurface, get_probabilities

#logging setting
logging.basicConfig(format='%(name)s  %(levelname)s  %(message)s', level=logging.INFO)

ads = ('H', )
cvgs = np.linspace(0, 1.0, 71)
colors = ['#000000', '#F08080', '#228B22', '#4169E1',
          '#EE7621', '#BA55D3', '#708090']
target_numbers = range(7)
probabilities_list = [[]]*7
for cvg in cvgs:
    surface = HexagonalSurface(10000, 10000)
    surface.initialize_surface(ads, (cvg, ), ads)
    probabilities = get_probabilities(surface, 'H', 'H', target_numbers)
    logging.info('coverage = %.3f, probabilities = %s', cvg, probabilities)
    # append colums(axis=1)
    probabilities = np.array(probabilities).reshape(-1, 1)
    probabilities_list = np.append(probabilities_list, probabilities, axis=1)

# plot
for idx, probabilities in enumerate(probabilities_list):
    plt.scatter(cvgs, probabilities, s=20, color=colors[idx],
                alpha=0.8, label=('nn = %d' % idx))
plt.xlabel(r'Coverage / ML')
plt.ylabel(r'Probability')
plt.grid(True)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.2)
plt.legend()

plt.show()
