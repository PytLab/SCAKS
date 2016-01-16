from operator import mul
import logging

import numpy as np
import matplotlib.pyplot as plt

from ..functions import convert_time


def get_square_circle_pts(nx, ny):
    '''
    Function to get positions of circles in grid,
    from left to right, from bottom to top.
    All number start from 1 instead of 0.

    Parameters:
    -----------
    nx: The number of circle along x axis, int

    ny: The number of circle along y axis, int

    Returns:
    --------
    map1d: a list of coordinates of grid,
           from left to right, from bottom to top,
           list of tuples of int.

    Example:
    --------
    >>> get_square_circle_pts(2, 2)
    >>> [(1, 1), (2, 1), (1, 2), (2, 2)]
    '''
    y, x = np.mgrid[1: nx+1, 1: ny+1]
    # use combination to get circle map
    map2d = [zip(a, b) for a, b in zip(x, y)]
    # reshape
    map1d = []
    for pts in map2d:
        map1d.extend(pts)

    return map1d


def get_square_lines_endpts(nx, ny):
    '''
    Function to get end points of square box line for plot use.

    Parameters:
    -----------
    nx: The number of circle along x axis, int
    ny: The number of circle along y axis, int

    Returns:
    --------
    endpts: a list of tuples, contains x and y endpts for a line.
            e.g. ((0.5, 2.5), (1.0, 1.0)) means x range is 0.5~2.5
            and y range is 1.0~1.0.

    Example:
    --------
    >>> get_square_lines_endpts(2,2)
    >>> [((0.5, 2.5), (1.0, 1.0)),
         ((0.5, 2.5), (2.0, 2.0)),
         ((1.0, 1.0), (0.5, 2.5)),
         ((2.0, 2.0), (0.5, 2.5))]

    '''
    # horizontal and vertical lines end points
    horizontal_endpts = [((0.5, ny+0.5), (float(y), float(y)))
                         for y in xrange(1, ny+1)]
    vertical_endpts = [((float(x), (float(x))), (0.5, nx+0.5))
                       for x in xrange(1, nx+1)]

    endpts = []
    endpts.extend(horizontal_endpts)
    endpts.extend(vertical_endpts)

    for endpt in endpts:
        logging.debug(str(endpt))

    return endpts


def plot_grid(shape, types, possible_types, color_dict,
              time, step, grid_type='square', circle_attrs={}):
    '''
    Function to plot lattice grid.

    Parameters:
    -----------
    shape: The grid shape, tuple of two int.

    types: The site types at the lattice points as a list,
           list of strings.

    possible_types: A list of possible types, list of strings.

    color_dict: circle color for different types.
                e.g. {type: color_code}, dict.

    time: time for the configure, float.

    grid_type: type of grid, 'square' or 'hexagonal', str.

    circle_attrs: same as **kwargs in plt.patch()
                  (http://matplotlib.org/api/patches_api.html), dict.
    '''
    if grid_type == 'square':
        get_lines_endpts = get_square_lines_endpts
        get_circle_pts = get_square_circle_pts

    # check consistency of color_dict and types
    consistent = np.array([t in color_dict for t in possible_types]).all()
    if not consistent:
        raise ValueError('All possible type should be in color_dict.')

    # check total point number
    consistent = (len(types) == mul(*shape))
    if not consistent:
        raise ValueError('points number in grid is not equal to types number.')

    # get grid lines end points & circle positions
    endpts = get_lines_endpts(*shape)
    cirpts = get_circle_pts(*shape)

    # set containers
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')

    # loop through type and points to add circle to axes
    for t, pt in zip(types, cirpts):  # t, pt <-> type, circle_point
        color = color_dict[t]
        # get circle object
        cir = plt.Circle(pt, radius=0.4, fc=color, **circle_attrs)
        ax.add_patch(cir)

    # add lines
    for endpt in endpts:
        ax.plot(*endpt, color='#ABABAB', linestyle='dotted', linewidth=0.2)

    # set axes attrs
    ax.set_xlim(0.5, shape[0]+0.5)
    ax.set_ylim(0.5, shape[1]+0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    # attrs of axis
    for spine in ax.spines.values():
        spine.set_linestyle('dashed')
        spine.set_alpha(0.5)
        spine.set_color('#AAAAAA')
    # get proper time format
    if time < 1e-2:
        time = '%es' % time
    else:
        time = '%dh %dmin %.2fsec' % convert_time(time)
    ax.set_title('Step = %d  ( %s )' % (step, time))

    return fig

if __name__ == '__main__':

    shape = (10, 10)
    types = ["O", "O", "C", "O", "O", "C", "O", "C", "C", "C", "C", "C", "C", "O", "O", "C", "O", "C", "C",
             "C", "C", "C", "C", "C", "C", "C", "C", "C", "O", "O", "C", "C", "C", "C", "C", "C", "C", "C",
             "C", "C", "C", "O", "O", "C", "C", "O", "O", "C", "O", "O", "C", "V", "C", "C", "C", "C", "C",
             "C", "C", "C", "C", "C", "C", "C", "C", "O", "O", "V", "C", "C", "O", "C", "C", "C", "C", "C",
             "C", "C", "C", "O", "C", "C", "C", "C", "O", "O", "C", "C", "O", "C", "C", "C", "C", "O", "O",
             "C", "C", "C", "O", "V"]
    possible_types = ('V', 'O', 'C')
    color_dict = dict(
        V='#FFFFFF',
        O='#FF6347',
        C='#607B8B',
        )
    circle_attrs = dict(
        alpha=0.7,
        antialiased=True,
        edgecolor='#FFFFFF',
        )

    fig = plot_grid(shape, types, possible_types, color_dict,
                    time=0.01, step=2, circle_attrs=circle_attrs)

    plt.show()
