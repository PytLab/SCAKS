import logging
import random


#logging setting
logging.basicConfig(format='%(levelname)s  %(message)s', level=logging.INFO)


class SurfaceSpaceError(Exception):
    "Exception raised for errors about surface space."
    pass


class LatticeSurface(object):
    def __init__(self, m, n):
        '''
        Basic solid lattice surface class,
        which could be inherited by other surface classes.
        '''
        self.shape = (m, n)

        self.grid = []
        for i in xrange(m):
            row = [0]*n
            self.grid.append(row)

        self.register_counter = {}

    def __str__(self):
        "Representation for grid on screen."
        s = ''
        m, n = self.shape

        for i in xrange(m):
            for j in xrange(n):
                if not self.grid[i][j]:
                    s += ('%-4s' % '.')
                else:
                    s += ('%-4s' % str(self.grid[i][j].symbol))
            s += '\n'

        return s

    def register(self, adsorbate):
        "Regist adsorbate to surface."
        x, y = adsorbate.x, adsorbate.y
        if not self.grid[x][y]:
            self.grid[x][y] = adsorbate
        else:
            raise SurfaceSpaceError('(%d, %d) is occupied by %s already' %
                                    (x, y, adsorbate.name))
        # add register counter
        if adsorbate.name in self.register_counter:
            self.register_counter[adsorbate.name] += 1
        else:
            self.register_counter.setdefault(adsorbate.name, 1)

        return


class HexagonalSurface(LatticeSurface):
    def __init__(self, m, n):
        LatticeSurface.__init__(self, m, n)

    def __str__(self):
        s = ''
        m, n = self.shape

        for i in xrange(m):
            if i % 2 != 0:  # translation
                s += ' '*2
            for j in xrange(n):
                if not self.grid[i][j]:
                    s += ('%-4s' % '.')
                else:
                    s += ('%-4s' % str(self.grid[i][j].symbol))
            s += '\n'

        return s

    def initialize_surface(self, ads, cvgs, symbols):
        '''
        Initialize surface randomly according to adsorbate coverage.

        Parameters:
        -----------
        ads: tuple of strings, names of adsorbates.
        cvgs: tuple of floats, coverages of adsorbates.
        symbols: tuple of strings, symbols of adsorbates.
        '''
        if len(ads) != len(cvgs):
            raise ValueError('shapes of ads and cvgs are not matched.')
        if not symbols:
            symbols = ads

        # get coverage ranges, [0.1, 0.2, 0.4, 0.3] -> [0.1, 0.3, 0.7, 1.0]
        l = len(cvgs)
        ranges = [sum(cvgs[: idx+1]) for idx in xrange(l)]
        ranges = [0.0] + ranges

        def find_adsorbate(ranges, p):
            for idx, r in enumerate(ranges):
                if p <= r:
                    return idx-1
            return -1

        m, n = self.shape
        # loop to fill surface with adsorbates
        for i in xrange(m):
            for j in xrange(n):
                p = random.random()
                # decide which adsorbate
                idx = find_adsorbate(ranges, p)
                if idx < 0:  # no adsorbate
                    continue
                adsorbate = Adsorbate(ads[idx], self, i, j, symbols[idx])
                self.register(adsorbate)
        return

    def get_NN_indices(self, x, y):
        '''
        Get the nearest neighbors indices tuple.
        '''
        # debug
        logging.debug('x, y = (%d, %d)', x, y)
        m, n = self.shape
        # NN positions order:
        #       1   2
        #     0   x   3
        #       5   4
        if x % 2 != 0:
            offsets = ((0, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0))
        else:
            offsets = ((0, -1), (-1, -1), (-1, 0), (0, 1), (1, 0), (1, -1))
        # debug
        logging.debug('offset = %s', str(offsets))

        # get nn indices
        def get_idx(offset):
            idx = [i + j for i, j in zip((x, y), offset)]
            # debug
            logging.debug('idx before pbc: %s', str(idx))
            # periodic boundary condition
            # row
            if idx[0] < 0:
                idx[0] = m - 1
            if idx[0] > m - 1:
                idx[0] = 0
            # column
            if idx[1] < 0:
                idx[1] = n - 1
            if idx[1] > n - 1:
                idx[1] = 0
            # debug
            logging.debug('idx after pbc: %s', str(idx))
            return tuple(idx)

        indices = [get_idx(offset) for offset in offsets]

        return tuple(indices)


class Adsorbate(object):
    def __init__(self, name, surface, x, y, symbol):
        '''
        Basic class for adsorbates adsorbed on surfaces.
        '''
        self.surface = surface
        self.name = name
        self.x, self.y = x, y
        self.symbol = symbol

    def __str__(self):
        return self.name
