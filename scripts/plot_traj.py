import os
import sys

from plot_grids import *

try:
    from PIL import Image
    PIL_installed = True
except ImportError:
    logging.info('PIL not installed, no animated gif would be created.')
    PIL_installed = False

import images2gif

# custom parameters
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

# check traj file existance
if len(sys.argv) > 1:
        filename = str(sys.argv[1])
else:
    if os.path.exists('./auto_trajectory.py'):
        filename = 'auto_trajectory.py'
    else:
        logging.error('No trajectory file supplied.')

# load traj info
globs = {}
locs = {}
execfile(filename, globs, locs)

# get grid shape
max_site = np.max(np.array(locs['sites']), axis=0)
shape = tuple([int(i) + 1 for i in max_site[:2]])
logging.debug('shape = %s', str(shape))

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
                    circle_attrs=circle_attrs)
    if not os.path.exists(path):
        os.mkdir(path)
    fname = path + str(step)+'.png'

    logging.info('creating %s ...', fname)
    fig.savefig(fname)
    logging.info('Ok.')

    # gather pic to generate animated gif
    if PIL_installed:
        im = Image.open(fname)
        images.append(im)

# make animated gif
if PIL_installed:
    gif_name = path + 'traj_movie.gif'
    logging.info('creating %s ...', gif_name)
    images2gif.writeGif(gif_name, images, duration=0.2)
    logging.info('Ok.')
