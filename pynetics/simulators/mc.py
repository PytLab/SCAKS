import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import *


def count_neighbors(surface, x, y, name):
    neighbor_indices = surface.get_NN_indices(x, y)
    counter = 0
    for nei_idx in neighbor_indices:
        i, j = nei_idx
        if surface.grid[i][j] and surface.grid[i][j].name == name:
            counter += 1

    return counter


def get_probabilities(surface, center_name, nei_name, target_numbers):
    '''
    Go through grid, to get probabilities of
    center_name has target_numbers nei_name in nearest neighbor.
    return a list of probabilities.
    '''
    length = len(target_numbers)
    if center_name not in surface.register_counter:
        return [0.0]*length

    m, n = surface.shape
    nn_counters = [0]*length
    for i in xrange(m):
        for j in xrange(n):
            adsorbate = surface.grid[i][j]
            if adsorbate and adsorbate.name == center_name:
                nei_number = count_neighbors(surface, i, j, nei_name)
                nn_counters[nei_number] += 1
    total_number = surface.register_counter[center_name]
    probabilities = [float(nn_counter)/total_number
                     for nn_counter in nn_counters]
    return probabilities


if __name__ == '__main__':
    ads = ('H', )
    cvgs = np.linspace(0, 1.0, 71)
    colors = ['#000000', '#F08080', '#228B22', '#4169E1', '#EE7621', '#BA55D3', '#708090']
    target_numbers = range(7)

    probabilities_list = [[]]*7
    for cvg in cvgs:
        surface = HexagonalSurface(100, 100)
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
    plt.xlabel('Coverage / ML')
    plt.ylabel('Probability')
    plt.grid(True)
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.2)
    plt.legend()

    plt.show()
