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
            filename = 'auto_trajectory.py'
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
