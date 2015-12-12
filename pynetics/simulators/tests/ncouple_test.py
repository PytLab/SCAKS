from mc_simulator import HexagonalSurface
import operator

size = (1000, 10000)

for i in xrange(50):
    surface = HexagonalSurface(*size)
    surface.initialize_surface(('H', 'CO'), (0.3, 0.7), ('H', 'CO'))
    ncouple = surface.count_couples(('H', 'CO'))
    print float(ncouple)/(operator.mul(*size))
