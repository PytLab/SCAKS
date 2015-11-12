import numpy as np
import matplotlib.pyplot as plt

from mc_simulator import C

cvgs = np.linspace(0, 1.0, 101)
colors = ['#000000', '#F08080', '#228B22', '#4169E1',
          '#EE7621', '#BA55D3', '#708090']

for i in xrange(7):
    p = [C(i, 6)*(theta**i)*((1 - theta)**(6 - i)) for theta in cvgs]
    plt.scatter(cvgs, p, s=20, color=colors[i], alpha=0.7, label=('nn = %d' % i))
    print C(i, 6)

plt.xlabel(r'Coverage / ML')
plt.ylabel(r'Probability')
plt.grid(True)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.2)
plt.legend()

plt.show()
