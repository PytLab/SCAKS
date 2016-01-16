'''
    Iterators for system of nonlinear equations roots finding.
'''

class RootfindingIterator(object):
    def __init__(self, f, x0, **kwargs):
        pass

    def __iter__(self):
        pass


class NewtonRoot(object):
    """
    Hacked from MDNewton in mpmath/calculus/optimization.py in order
    to allow for constraints on the solution.
    Find the root of a vector function numerically using Newton's method.

    f is a vector function representing a nonlinear equation system.

    x0 is the starting point close to the root.

    J is a function returning the Jacobian matrix for a point.

    constraint is function to limit x.
    """
    def __init__(self, f, J, x0, constraint, norm,
                 mpfloat, matrix, Axb_solver, **kwargs):
        self.f, self.x0, self.J = f, x0, J
        #below are all function objects
        self._norm = norm
        self._mpf = mpfloat
        self._matrix = matrix
        self._Axb_solver = Axb_solver

#        if 'constraint' in kwargs:
#            self.constraint = kwargs['constraint']
#        else:
#            def constraint(x):
#                return x
#            self.constraint = constraint
        self.real_constraint = constraint

        def quasi_constraint(x):
            return x
        self.quasi_constraint = quasi_constraint

        if 'dtheta_dt_expressions' in kwargs:
            self.dtheta_dt_expressions = kwargs['dtheta_dt_expressions']

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
        J = self.J
        x0 = self.quasi_constraint(self.x0)
        norm = self._norm
        fx = self._matrix(f(x0))
        fxnorm = norm(fx)
        cancel = False
        #x0 = mp.matrix(x0)
        while not cancel:
            iter_counter += 1
            if iter_counter <= 5:
                self.constraint = self.quasi_constraint
            else:
                self.constraint = self.real_constraint
            #get direction of descent
            fx = self._matrix(f(x0))
            fxn = -fx
            Jx = J(self.dtheta_dt_expressions, x0)
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
