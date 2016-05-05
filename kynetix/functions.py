# -*- coding: utf-8 -*-
#↓↓↓static utility functions↓↓↓
import copy


def string2symbols(s):  # HACKED# copy-pasted from ase\atoms.py
    """Convert string to list of chemical symbols."""
    n = len(s)

    if n == 0:
        return []

    c = s[0]

    if c.isdigit():
        i = 1
        while i < n and s[i].isdigit():
            i += 1
        return int(s[:i]) * string2symbols(s[i:])

    if c == '(':
        p = 0
        for i, c in enumerate(s):
            if c == '(':
                p += 1
            elif c == ')':
                p -= 1
                if p == 0:
                    break
        j = i + 1
        while j < n and s[j].isdigit():
            j += 1
        if j > i + 1:
            m = int(s[i + 1:j])
        else:
            m = 1
        return m * string2symbols(s[1:i]) + string2symbols(s[j:])

    if c.isupper():
        i = 1
        if 1 < n and s[1].islower():
            i += 1
        j = i
        while j < n and s[j].isdigit():
            j += 1
        if j > i:
            m = int(s[i:j])
        else:
            m = 1
        return m * [s[:i]] + string2symbols(s[j:])
    else:
        raise ValueError


def numerical_jacobian(f, x, matrix, num_repr='gmpy', h=1e-10, direction='right'):
    """
    Calculate the Jacobian matrix of a function at the point x0.

    This is the first derivative of a vectorial function:

        f : R^m -> R^n with m >= n

    Direction is the direction of h (perturbation).

    Modified from CatMAP.
    """
    def vec2list(col_vector):
        "convert column vector of mpmath or numpy to python list."
        if num_repr == 'gmpy':
            return tuple(col_vector.reshape(1, -1).tolist()[0])
        elif num_repr == 'mpmath':
            return col_vector

    #x = matrix(x)
    fx = matrix(f(x))
    m = len(fx)
    n = len(x)
    J = matrix(m, n)
    for j in xrange(n):
        xj = copy.copy(x)
        delta = abs(h*xj[j])
        delta = max(delta, h)
        #using delta proportional to xj is more stable
        #for very small numbers.
        if direction == 'left':
            xj[j] -= delta
            Jj = (matrix(f(xj)) - fx)/(-delta)
            Jj = vec2list(Jj)
        if direction == 'right':
            xj[j] += delta
            Jj = (matrix(f(xj)) - fx)/(delta)
            Jj = vec2list(Jj)
        for i in xrange(m):
            J[i, j] = Jj[i]
    return J


def get_list_string(var_name, list_obj):
    '''
    Function to get string format of a list object.

    Parameters:
    -----------
    var_name: name of variable, str.

    list_obj: object of variable.

    Example:
    --------
    >>> a = [(1,2,3), (2,3,4)]
    >>> get_list_string('a', a)
    >>> 'a = [\n    (1, 2, 3),\n    (2, 3, 4),\n]\n\n'

    Called:
    -------
    kmc_plugins.CoveragesAnalysis()

    '''
    begin = var_name + ' = ['
    indent = ' '*4
    data = ''
    for idx, elem in enumerate(list_obj):
        # if item is iterable, one a line
        if hasattr(elem, '__iter__'):
            data += ('\n' + indent + str(elem) + ',')
            continue
        # 5 items a line by default
        if idx % 5 == 0:
            data += ('\n' + indent)

        # add single quotes for string
        if isinstance(elem, str):
            data += ("'" + elem + "', ")
        else:
            data += (str(elem) + ', ')
    end = '\n]\n\n'

    content = begin + data + end

    return content


def convert_time(sec):
    '''
    Convert format of time from seconds to *h *min *sec.
    '''
    hours = int(sec/(3600.0))
    minutes = int((sec - hours*3600)/60.0)
    seconds = sec - hours*3600 - minutes*60

    return hours, minutes, seconds  # int, int, float

