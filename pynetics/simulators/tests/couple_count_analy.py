from math import sqrt

import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import C


def couple_probability(thetaA, thetaB, thetaF):
    def f1(theta):
        return 3*theta*(1 - theta)**2

    def f2(theta):
        return 3*theta**2*(1 - theta)

    def f3(theta):
        return theta**3

    def f0(theta):
        return (1 - theta)**3

    functions = [f0, f1, f2, f3]

    p_tot = 2*thetaA*thetaB
    p_correct = 0.0
    for i, alpha in enumerate(functions):
        for j, beta in enumerate(functions):
            p_correct += (i*j+1)*(alpha(thetaA)*beta(thetaB) *
                          (thetaF**2/((i+1)*(j+1)) + thetaA**2/((j+1)*(i+3)) +
                           thetaB**2/((i+1)*(j+3)) + 2*thetaA*thetaB/((i+2)*(j+2))))

    return p_tot*p_correct


colors = ['#000000', '#F08080', '#228B22', '#4169E1',
          '#EE7621', '#BA55D3', '#708090', '#FFFF00', '#EE00EE']
co_cvg_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
line_list = []
h_cvgs_list = []

for co_cvg in co_cvg_list:
    h_cvgs = np.linspace(0, 1 - co_cvg, 71)
    h_cvgs_list.append(h_cvgs)
    probabilities = []
    for h_cvg in h_cvgs:
#        p = h_cvg*co_cvg
        p = couple_probability(co_cvg, h_cvg, 0.0)
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
