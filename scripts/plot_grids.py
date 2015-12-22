from operator import mul
import logging

import numpy as np
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO)


def get_square_circle_pts(nx, ny):
    '''
    Function to get positions of circles in grid,
    from left to right, from bottom to top.
    All number start from 1 instead of 0.

    Parameters:
    -----------
    :param nx: The number of circle along x axis.
    :type nx: int

    :param ny: The number of circle along y axis.
    :type ny: int
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
    :param nx: The number of circle along x axis.
    :type nx: int
    :param ny: The number of circle along y axis.
    :type ny: int
    '''
    # horizontal and vertical lines end points
    horizontal_endpts = [((0.5, ny+0.5), (float(y), float(y)))
                         for y in xrange(1, ny+1)]
    vertical_endpts = [((float(x), (float(x))), (0.5, nx+0.5))
                       for x in xrange(1, nx+1)]

    endpts = []
    endpts.extend(horizontal_endpts)
    endpts.extend(vertical_endpts)

    # debug
    logging.debug('line end points = ')
    for endpt in endpts:
        logging.debug(str(endpt))

    return endpts


def plot_grid(shape, types, possible_types, color_dict,
              time, step, grid_type='square', circle_attrs={}):
    '''
    Function to plot lattice grid.

    Parameters:
    -----------
    :param shape: The grid shape.
    :type shape: tuple of two int.

    :param types: The site types at the lattice points as a list.
    :type types: list of strings.

    :param possible_types: A list of possible types.
    :type possible_types: list of strings.

    :param color_dict: circle color for different types.
                       e.g. {type: color_code}
    :type color_dict: dict.

    :param time: time for the configure.
    :type time: float.

    :param grid_type: type of grid, 'square' or 'hexagonal'
    :type grid_type: str.

    :param **circle_attrs: same as **kwargs in plt.patch()
                           (http://matplotlib.org/api/patches_api.html)
    :type **circle_attrs: dict.
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
        ax.plot(*endpt, color='#ABABAB', linestyle='dashed', linewidth=1.2)

    # set axes attrs
    ax.set_xlim(0.5, shape[0]+0.5)
    ax.set_ylim(0.5, shape[1]+0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Step = %d   Time = %.2e' % (step, time))

    return fig

if __name__ == '__main__':

    shape = (10, 10)
    types = ["O","O","C","O","O","C","O","C","C","C","C","C","C","O","O","C","O","C","C",
             "C","C","C","C","C","C","C","C","C","O","O","C","C","C","C","C","C","C","C",
             "C","C","C","O","O","C","C","O","O","C","O","O","C","V","C","C","C","C","C",
             "C","C","C","C","C","C","C","C","O","O","V","C","C","O","C","C","C","C","C",
             "C","C","C","O","C","C","C","C","O","O","C","C","O","C","C","C","C","O","O",
             "C","C","C","O","V"]
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
