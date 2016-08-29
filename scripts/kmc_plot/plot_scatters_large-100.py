import os
import logging
logging.basicConfig(level=logging.INFO)
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
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

from kynetix.utilities.format_utilities import convert_time
from kynetix.plotters import images2gif


def plot_scatters(types,
                  shape,
                  coordinates,
                  possible_types,
                  color_dict,
                  time, step,
                  circle_attrs=None):
    """
    Function to plot lattice grid.

    Parameters:
    -----------
    types: The site types at the lattice points as a list,
           list of strings.

    possible_types: A list of possible types, list of strings.

    color_dict: circle color for different types.
                e.g. {type: color_code}, dict.

    time: time for the configure, float.

    circle_attrs: same as **kwargs in plt.patch()
                  (http://matplotlib.org/api/patches_api.html), dict.
    """
    # Check consistency of color_dict and types.
    consistent = np.array([t in color_dict for t in possible_types]).all()
    if not consistent:
        raise ValueError('All possible type should be in color_dict.')

    # Set containers.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')

    # Get scatter attrs
    alpha = circle_attrs['alpha'] if 'alpha' in circle_attrs else 0.7
    if "area" not in circle_attrs:
        area = 60.0**2/len(types)*20
    else:
        area = circle_attrs['area']

    # Classify points.
    scatter_dict = {}  # {adsorbate string: list of points}
    for t in possible_types:
        scatter_dict.setdefault(t, [])

    for t, coord in zip(types, coordinates):
        scatter_dict[t].append(coord[: 2])

    for t, pts in scatter_dict.iteritems():
        if not pts:
            continue
        x, y = zip(*pts)
        color = color_dict[t]
        alpha = circle_attrs['alpha'] if 'alpha' in circle_attrs else 0.7
        edgecolor = circle_attrs['edgecolor'] if 'edgecolor' in circle_attrs else color
        marker = circle_attrs['marker'] if 'marker' in circle_attrs else 's'
        ax.scatter(x, y, s=area, c=color, alpha=alpha, edgecolor=edgecolor, marker=marker)

    # Set axes attrs.
    ax.set_xlim(-0.5, shape[0])
    ax.set_ylim(-0.5, shape[-1])
    ax.set_xticks([])
    ax.set_yticks([])

    # Attrs of axis.
    for spine in ax.spines.values():
        spine.set_linestyle('dashed')
        spine.set_alpha(0.5)
        spine.set_color('#AAAAAA')

    # Get proper time format
    if time < 1e-2:
        time = '%es' % time
    else:
        time = '%dh %dmin %.2fsec' % convert_time(time)
    ax.set_title('Step = %d  ( %s )' % (step, time))

    return fig

if __name__ == '__main__':

    # Set argument parser.
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gif",
                        help="create gif animated picture.",
                        action="store_true")
    args = parser.parse_args()

    shape = (60, 60)

    possible_types = ("O_u", "O_d", "O_l", "O_r", "V", "O_s", "C")

    color_dict = dict(
        V='#FFFFFF',
        O_s='#EE0000',
        O_u='#FF6347',
        O_d='#FF6347',
        O_l='#FF6347',
        O_r='#FF6347',
        C='#363636',  # '#607B8B'
        )
    circle_attrs = dict(
        area=2.08,
        alpha=0.7,
        #antialiased=True,
        )

    # Locate trajectory file.
    filename = 'auto_lattice_trajectory.py'
    if not os.path.exists(filename):
        raise IOError('No trajectory file found.')

    # Read data from file.
    globs = {}
    locs = {}
    execfile(filename, globs, locs)

    steps = locs["steps"]
    times = locs["times"]
    sites = locs["sites"]
    types = locs["types"]
    images = []
    path = "./trajplots/"

    for tp, step, simu_time in zip(types, steps, times):
        fig = plot_scatters(tp, shape, sites, possible_types, color_dict,
                            time=simu_time, step=step, circle_attrs=circle_attrs)
        if not os.path.exists(path):
                os.mkdir(path)
        fname = path + str(step) + '.png'

        logging.info("creating {} ...".format(fname))
        fig.savefig(fname)

        # Clear figure object after saving.
        plt.close(fig)
        logging.info("Ok.")

        if PIL_installed:
            im = Image.open(fname)
            images.append(im)

    # Make gif.
    if PIL_installed and args.gif:
        gif_name = path + 'traj_movie.gif'
        logging.info('creating %s ...', gif_name)
        images2gif.writeGif(gif_name, images, duration=0.3)
        logging.info('Ok.')
