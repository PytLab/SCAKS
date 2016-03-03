import logging
import os

import numpy as np
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
from .grid_plot import *
from .scatter_plot import *


class KMCPlotter(PlotterBase):
    def __init__(self, owner):
        '''
        Class for kMC map plot.
        '''
        PlotterBase.__init__(self, owner)
        # set logger
        self.logger = logging.getLogger('model.plotters.KMCPlotter')

    def get_parameters(func):
        '''
        Decorator to extract attrs values from model,
        and set them as parameters for kmc plotting.
        '''
        def wrapped_func(self, types, time, step, circle_attrs={}):
            # get essential parameters
            shape = self._owner.grid_shape
            adsorbate_names = [ads.split('_')[0] for ads in self._owner.adsorbate_names]
            possible_types = adsorbate_names + ['Vac']
            color_dict = self._owner.color_dict
            grid_type = self._owner.grid_type
            circle_attrs = self._owner.circle_attrs

            fig = func(self, types, time, step, circle_attrs={})

            return fig

        return wrapped_func

    @get_parameters
    def plot_grid(self, types, time, step, circle_attrs={}):
        '''
        Function to plot lattice grid of a specified configure.

        Returns:
        --------
        fig: plt.figure object contains grid plot.

        '''
        # parameters are get by decorator @get_parameters

        # plot grid
        fig = plot_grid(shape, types, possible_types, color_dict,
                        time, step, grid_type, circle_attrs)

        return fig

    @get_parameters
    def plot_scatter(self, types, time, step, circle_attrs={}):
        '''
        Function to plot lattice grid of a specified configure using scatter plot.

        Returns:
        --------
        fig: plt.figure object contains grid plot.

        '''
        # parameters are get by decorator @get_parameters

        # plot grid
        fig = plot_scatter(shape, types, possible_types, color_dict,
                           time, step, grid_type, circle_attrs)

        return fig

    def plot_traj(self, mode='grid', gif=True, filename=None):
        '''
        Function plot trajectory pictures and create animated gif.

        Parameters:
        -----------
        mode: grid plotting mode, | 'grid' | 'scatter' |, defaults to be 'grid', str.

        gif: generate GIF animation or not, defaults to be True, bool.

        filename: input trajectory filename, use 'auto_trajectory.py' to be default, str.

        Example:
        >>> m.plotter.plot_traj(mode='scatter', gif=False, filename='traj.py')

        '''
        # get essential parameters
        shape = self._owner.grid_shape
        adsorbate_names = [ads.split('_')[0] for ads in self._owner.adsorbate_names]
        possible_types = adsorbate_names + ['Vac']
        color_dict = self._owner.color_dict
        grid_type = self._owner.grid_type
        circle_attrs = self._owner.circle_attrs

        # locate trajectory file
        if not filename:
            filename = 'auto_trajectory.py'
        if not os.path.exists(filename):
            self.logger.error('No trajectory file found.')
            raise FilesError('No trajectory file found.')

        # get plot function
        if mode == 'grid':
            plot_func = plot_grid
        elif mode == "scatter":
            plot_func = plot_scatters
        else:
            raise ParameterError('mode must be \'grid\' or \'scatter\'.')

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
            fig = plot_func(shape, tp, possible_types, color_dict,
                            time=simu_time, step=step,
                            grid_type=grid_type,
                            circle_attrs=circle_attrs)
            if not os.path.exists(path):
                os.mkdir(path)
            fname = path + str(step) + '.png'

            self.logger.info('creating %s ...', fname)
            fig.savefig(fname)
            self.logger.info('Ok.')

            # gather pic to generate animated gif
            if PIL_installed:
                im = Image.open(fname)
                images.append(im)

        # make animated gif
        if PIL_installed and gif:
            gif_name = path + 'traj_movie.gif'
            self.logger.info('creating %s ...', gif_name)
            images2gif.writeGif(gif_name, images, duration=0.3)
            self.logger.info('Ok.')
