'''
    Iterators for system of nonlinear equations roots finding.
'''
from ..errors.error import *


class RootfindingIterator(object):
    def __init__(self, f, x0, **kwargs):
        # set attributes
        self.f, self.x0 = f, x0        
        for key in kwargs:
            setattr(self, '_' + key, kwargs[key])

    def __iter__(self):
        pass


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
        for param in ['J', 'constraint', 'norm', 'mpfloat', 'matrix',
                      'Axb_solver', 'dtheta_dt_expressions']:
            if not hasattr(self, '_' + param):
                msg = 'parameter \'%s\' must be supplied.' % param
                raise ParameterError(msg)

        # set constraint function
        self.real_constraint = self._constraint

        def pseudo_constraint(x):
            return x
        self.pseudo_constraint = pseudo_constraint

    def __iter__(self):

        def vec2tup(col_vector):
            "convert column vector of mpmath or numpy to python tuple."
            if self._mpf == gmpy2.mpfr:  # gmpy
                return tuple(col_vector.reshape(1, -1).tolist()[0])
#            elif self._mpf == mp.mpf:  # mpmath
#                return tuple(col_vector)
#            elif self._mpf == sym.mpmath.mpf:  # sympy
#                return tuple(col_vector)
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
            Jx = J(self._dtheta_dt_expressions, x0)
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
                print "Solver: Found stationary point."
                cancel = True
            fx = self._matrix(f(x1))
            x0, fxnorm = x1, norm(fx)

            yield x0, fxnorm, fx
