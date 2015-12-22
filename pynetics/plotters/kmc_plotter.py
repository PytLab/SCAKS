import logging
from operator import mul
import os

import numpy as np
import matplotlib.pyplot as plt
try:
    from PIL import Image
    PIL_installed = True
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                   !!!"
    print "!!!            WARNING: PIL is not installed          !!!"
    print "!!!          No animated gif would be created.        !!!"
    print "!!!                                                   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    PIL_installed = False

import images2gif
from plotter_base import PlotterBase
from ..errors.error import *


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
    ax.set_title('Step = %d   Time = %fs' % (step, time))

    return fig


class KMCPlotter(PlotterBase):
    def __init__(self, owner):
        '''
        Class for kMC map plot.
        '''
        PlotterBase.__init__(self, owner)
        # set logger
        self.logger = logging.getLogger('model.plotters.KMCPlotter')

    def plot_grid(self, types, time, step, circle_attrs={}):
        '''
        Function to plot lattice grid of a specified configure.

        Returns:
        --------
        fig: plt.figure object contains grid plot.

        '''
        # get essential parameters
        shape = self._owner.grid_shape
        adsorbate_names = [ads.split('_')[0] for ads in self._owner.adsorbate_names]
        possible_types = adsorbate_names + ['Vac']
        color_dict = self._owner.color_dict
        grid_type = self._owner.grid_type
        circle_attrs = self._owner.circle_attrs

        # plot grid
        fig = plot_grid(shape, types, possible_types, color_dict,
                        time, step, grid_type, circle_attrs)

        return fig

    def plot_traj(self, filename=None):
        '''
        Function plot trajectory pictures and create animated gif.
        '''
        # get essential parameters
        adsorbate_names = [ads.split('_')[0] for ads in self._owner.adsorbate_names]
        possible_types = adsorbate_names + ['Vac']
        color_dict = self._owner.color_dict
        grid_type = self._owner.grid_type
        circle_attrs = self._owner.circle_attrs

        # locate trajectory file
        if not filename:
            filename = 'trajectory.py'
        if not os.path.exists(filename):
            self.logger.error('No trajectory file found.')
            raise FilesError('No trajectory file found.')

        # load traj info
        globs = {}
        locs = {}
        execfile(filename, globs, locs)

        # get grid shape
        max_site = np.max(np.array(locs['sites']), axis=0)
        shape = tuple([int(i) + 1 for i in max_site[:2]])
        self.logger.info('shape = %s', str(shape))

        # get steps and times
        steps = locs['steps']
        times = locs['times']

        # go through types to plot grids
        types = locs['types']
        images = []
        path = './trajplots/'
        for tp, step, simu_time in zip(types, steps, times):
            fig = plot_grid(shape, tp, possible_types, color_dict,
                            time=simu_time, step=step,
                            grid_type=grid_type,
                            circle_attrs=circle_attrs)
            if not os.path.exists(path):
                os.mkdir(path)
            fname = path + str(step)+'.png'

            self.logger.info('creating %s ...', fname)
            fig.savefig(fname)
            self.logger.info('Ok.')

            # gather pic to generate animated gif
            if PIL_installed:
                im = Image.open(fname)
                images.append(im)

        # make animated gif
        if PIL_installed:
            gif_name = path + 'traj_movie.gif'
            self.logger.info('creating %s ...', gif_name)
            images2gif.writeGif(gif_name, images, duration=0.2)
            self.logger.info('Ok.')
