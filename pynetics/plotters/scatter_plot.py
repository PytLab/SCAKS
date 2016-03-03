from operator import mul

import numpy as np
import matplotlib.pyplot as plt

from ..functions import convert_time


def get_square_scatter_pts(nx, ny):
    '''
    Function to get positions of scatters in grid,
    from left to right, from bottom to top.

    Parameters:
    -----------
    nx: The number of circle along x axis, int

    ny: The number of circle along y axis, int

    Returns:
    --------
    map1d: a list of coordinates of scatter,
           from left to right, from bottom to top,
           list of tuples of float.

    Example:
    --------
    >>> get_square_scatter_pts(2, 2)
    >>> [(0.0, 0.0), (0.0, 0.1), (0.1, 0.0), (0.1, 0.1)]
    '''
    # get possible x and y values
    xv = np.arange(0.0, 0.1*nx, 0.1)  # start from 0.0
    yv = np.arange(0.0, 0.1*ny, 0.1)

    # get all x values of points
    x = xv.repeat(ny, axis=0)

    # get all y valus of points
    yv = yv.reshape(-1, 1)
    y = yv.repeat(nx, axis=1)
    y = y.reshape(-1, order='F')

    return zip(x, y)


def plot_scatters(shape, types, possible_types, color_dict,
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
        get_scatter_pts = get_square_scatter_pts

    # check consistency of color_dict and types
    consistent = np.array([t in color_dict for t in possible_types]).all()
    if not consistent:
        raise ValueError('All possible type should be in color_dict.')

    # check total point number
    consistent = (len(types) == mul(*shape))
    if not consistent:
        raise ValueError('points number in grid is not equal to types number.')

    # get grid lines end points & circle positions
    scatter_pts = get_scatter_pts(*shape)

    # set containers
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')

    # classify points
    scatter_dict = {}  # {adsorbate string: list of points}
    for t in possible_types:
        scatter_dict.setdefault(t, [])

    for t, pt in zip(types, scatter_pts):  # t, pt <-> type, scatter_pt
        color = color_dict[t]
        scatter_dict[t].append(pt)

    # plot scatter points
    for t, pts in scatter_dict.iteritems():
        if not pts:
            continue
        x, y = zip(*pts)
        # get scatter attrs
        color = color_dict[t]
        alpha = circle_attrs['alpha'] if 'alpha' in circle_attrs else 0.7
        edgecolor = circle_attrs['edgecolor'] if 'edgecolor' in circle_attrs else color
        area = circle_attrs['area'] if 'area' in circle_attrs else 60.0**2/mul(*shape)*20
        marker = circle_attrs['marker'] if 'marker' in circle_attrs else 'o'
        ax.scatter(x, y, s=area, c=color, alpha=alpha, edgecolor=edgecolor, marker=marker)

    # set axes attrs
    ax.set_xlim(-0.1, shape[0]*0.1)
    ax.set_ylim(-0.1, shape[-1]*0.1)
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
    ax.grid(True)

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
        marker='s',
        area=20,
        alpha=0.7,
        antialiased=True,
        edgecolor='#FFFFFF',
        )

    fig = plot_scatters(shape, types, possible_types, color_dict,
                        time=0.01, step=2, circle_attrs=circle_attrs)

    plt.show()
