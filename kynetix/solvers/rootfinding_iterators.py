'''
    Iterators for system of nonlinear equations roots finding.
'''

import logging

import mpmath as mp
from scipy.optimize import golden
import gmpy2
import sympy as sym

from kynetix.errors.error import *


class RootfindingIterator(object):
    '''
    Class base for rootfinding iterator sub-classes to inherit from.
    '''
    def __init__(self, f, x0, **kwargs):
        # set attributes
        self.f, self.x0 = f, x0
        for key in kwargs:
            setattr(self, '_' + key, kwargs[key])

    def __iter__(self):
        '''
        Returns:
        -------
        x0: current x vector, tuple of float.

        fxnorm: norm of f(x), float

        fx: value of f(x), sequence of float
        '''
        pass
        return (x0, fxnorm, fx)


class ConstrainedNewton(RootfindingIterator):
    """
    Hacked from MDNewton in mpmath/calculus/optimization.py in order
    to allow for constraints on the solution.
    Find the root of a vector function numerically using Newton's method.

    Parameters:
    -----------
    f is a vector function representing a nonlinear equation system.

    x0 is the starting point close to the root.

    **kwargs MUST contains: **
    J: a function returning the Jacobian matrix for a point.

    constraint: a coverages tuple constraint function to set limit to x.

    norm: a function to get a norm.

    mpfloat: high-precision float type, e.g. mpmath.mpfloat

    matrix: matrix type, e.g. mpmath.matrix

    Axb_solver: a function to solve system of linear equations by solving Ax=b.

    """
    def __init__(self, f, x0, **kwargs):
        RootfindingIterator.__init__(self, f, x0, **kwargs)

        # check essential parameters reading
        essential_params = ('J', 'constraint', 'norm', 'mpfloat',
                            'matrix', 'Axb_solver')
        for param in essential_params:
            if not hasattr(self, '_' + param):
                msg = "parameter '{}' must be supplied.".format(param)
                raise ParameterError(msg)

        # set constraint function
        self.real_constraint = self._constraint

        def pseudo_constraint(x):
            return x
        self.pseudo_constraint = pseudo_constraint

        # set logger
        self.logger = logging.getLogger('model.solvers.ConstrainedNewton')

    def __iter__(self):

        def vec2tup(col_vector):
            "convert column vector of mpmath or numpy to python tuple."
            if self._mpfloat == gmpy2.mpfr:  # gmpy
                return tuple(col_vector.reshape(1, -1).tolist()[0])
            else:
                return tuple(col_vector)

        iter_counter = 0
        f = self.f
        J = self._J
        x0 = self.pseudo_constraint(self.x0)
        norm = self._norm
        fx = self._matrix(f(x0))
        fxnorm = norm(fx)
        cancel = False
        #x0 = mp.matrix(x0)
        while not cancel:
            iter_counter += 1
            if iter_counter <= 5:
                self.constraint = self.pseudo_constraint
            else:
                self.constraint = self.real_constraint
            #get direction of descent
            fx = self._matrix(f(x0))
            fxn = -fx
            Jx = J(x0)
            try:
                s = self._Axb_solver(Jx, fxn)  # if use gmpy and numpy,
                #print s                       # lose precision here
            except ZeroDivisionError:
                #print 'ZeroDivisionError!'
                #cancel = True
                #break
                raise ValueError("ZeroDivisionError!")

            #use golden method to get optimal step size
            def fl(l):
                x1 = self._matrix(x0) + l*s
                fx = self._matrix(f(vec2tup(x1)))
                return norm(fx)
            l = golden(fl)
#            print l
#            l = mp.mpf('1.0')
            x1 = self._matrix(x0) + l*s  # matrix
            x1 = self.constraint(vec2tup(x1))
            if x1 == x0:
                self.logger.info("Solver: Found stationary point.")
                cancel = True
            fx = self._matrix(f(x1))
            x0, fxnorm = x1, norm(fx)

            yield (x0, fxnorm, fx)


class MDNewton(RootfindingIterator):
    """
    Find the root of a vector function numerically using Newton's method.

    f is a vector function representing a nonlinear equation system.

    x0 is the starting point close to the root.

    J is a function returning the Jacobian matrix for a point.

    Supports overdetermined systems.

    Use the 'norm' keyword to specify which norm to use. Defaults to max-norm.
    The function to calculate the Jacobian matrix can be given using the
    keyword 'J'. Otherwise it will be calculated numerically.

    Please note that this method converges only locally. Especially for high-
    dimensional systems it is not trivial to find a good starting point being
    close enough to the root.

    It is recommended to use a faster, low-precision solver from SciPy [1] or
    OpenOpt [2] to get an initial guess. Afterwards you can use this method for
    root-polishing to any precision.

    [1] http://scipy.org

    [2] http://openopt.org/Welcome
    """

    def __init__(self, f, x0, **kwargs):
        self.f = f
        if isinstance(x0, (tuple, list)):
            x0 = mp.matrix(x0)
        assert x0.cols == 1, 'need a vector'
        self.x0 = x0
        if 'J' in kwargs:
            self.J = kwargs['J']
        else:
            def J(*x):
                return mp.jacobian(f, x)
            self.J = J
        self.norm = mp.norm
        self.verbose = kwargs['verbose']

        # set logger
        self.logger = logging.getLogger('model.solvers.MDNewton')

    def __iter__(self):
        f = self.f
        x0 = self.x0
        norm = self.norm
        J = self.J
        fx = mp.matrix(f(x0))
        fxnorm = norm(fx)
        cancel = False
        while not cancel:
            # get direction of descent
            fxn = -fx
            Jx = J(x0)
            s = mp.lu_solve(Jx, fxn)
            if self.verbose:
                self.logger.debug('Jx = \n%s', str(Jx))
                self.logger.debug('s = \n%s', str(s))
            # damping step size TODO: better strategy (hard task)
            l = mp.mpf('1.0')
            x1 = x0 + s
            while True:
                if x1 == x0:
                    self.logger.info("Found stationary point.")
                    cancel = True
                    break
                fx = mp.matrix(f(x1))
                newnorm = norm(fx)
                if newnorm < fxnorm:
                    # new x accepted
                    fxnorm = newnorm
                    x0 = x1
                    break
                l /= 2
                x1 = x0 + l*s

            yield (tuple(x0), fxnorm, fx)

