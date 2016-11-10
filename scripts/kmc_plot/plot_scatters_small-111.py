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
                  markers,
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

    markers: A list of marker parameter for corresponding basis_site.

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
        area = circle_attrs["area"]

    # Classify points.
    scatter_dict = {}  # {adsorbate string: list of points}
    for t in possible_types:
        scatter_dict.setdefault(t, [])

    site_markers = ['o', 'v', 'v', 'v', 's', 's']
    for idx, (t, coord) in enumerate(zip(types, coordinates)):  # t, pt <-> type, scatter_pt
        x, y = coord[: 2]
        color = color_dict[t]
        #marker = markers[t]
        marker = site_markers[idx % 6]
        edgecolor = color
        ax.scatter(x, y, s=area, c=color, alpha=alpha, edgecolor=edgecolor, marker=marker)

    # Set axes attrs.
    ax.set_xlim(-0.5, 3./2*shape[0])
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

    shape = (10, 10)

    possible_types = ("O_u", "O_d", "O_l", "O_r", "O_ur", "O_dr", "O_dl", "O_ul", "V", "O_s", "C")

    markers = dict(V='o',
                   O_s='o',
                   O_u='^',
                   O_d='v',
                   O_l='<',
                   O_r='>',
                   O_ur='x',
                   O_dr='x',
                   O_ul='x',
                   O_dl='x',
                   C='o')

    color_dict = dict(V='#FFFFFF',
                      O_s='#FF6347',
                      O_u='#EE0000',
                      O_d='#EE0000',
                      O_l='#EE0000',
                      O_r='#EE0000',
                      O_ur='#EE0000',
                      O_dr='#EE0000',
                      O_ul='#EE0000',
                      O_dl='#EE0000',
                      C='#607B8B')

    circle_attrs = dict(alpha=0.7,
                        antialiased=True,
                        area=20,
                        edgecolor='#FFFFFF')

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

    # Get cartisian coordinates.
    cell_vectors = np.matrix([[1.0, 0.0, 0.0],
                              [0.5, 0.87, 0.0],
                              [0.0, 0.0, 1.0]])
    sites = np.matrix(sites)
    sites = (sites*cell_vectors).tolist()

    for tp, step, simu_time in zip(types, steps, times):
        fig = plot_scatters(tp, shape, sites, markers, possible_types, color_dict,
                            time=simu_time, step=step, circle_attrs=circle_attrs)
        if not os.path.exists(path):
            os.mkdir(path)
        fname = path + str(step) + '.png'

        if os.path.exists(fname):
            fname = "{}{}-redis.png".format(path, step)

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

