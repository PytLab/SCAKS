import logging
logging.basicConfig(level=logging.INFO)
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from kynetix.utilities.coordinate_utilities import CoordsGroup


def plot_elements(cell_vectors,
                  coords_group,
                  color_dict,
                  marker_dict,
                  title):
    """
    Function to plot CoordsGroup object.
    """
    # Get cartisan coordinates.
    dir_coords = np.matrix(coords_group.coordinates())
    cell_vectors = np.matrix(cell_vectors)
    cart_coords = np.array(dir_coords*cell_vectors)
    xs = cart_coords[:, 0].tolist()
    ys = cart_coords[:, 1].tolist()

    fig = plt.figure(frameon=False)
    ax = fig.add_subplot(111)
    ax.axis('off')
    ax.set_aspect("equal", adjustable="box")

    # Remove axis borders.
    #for item in [fig, ax]:
    #    item.patch.set_visible(False)

    # Plot scatters.
    types = coords_group.elements()
    for tp, x, y in zip(types, xs, ys):
        color = color_dict[tp]
        marker = marker_dict[tp]
        if color in ["#FFFFFF", "#FFF"]:
            edgecolor = "#5E5E5E"
        else:
            edgecolor = color
        ax.scatter(x, y, s=500, c=color, alpha=0.7,
                   edgecolor=edgecolor, marker=marker)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(min(xs)-0.5, max(xs)+0.5)
    ax.set_ylim(min(ys)-0.5, max(ys)+0.5)
    ax.set_title(title)

    return fig


def plot_process(process_dict,
                 basis_coords,
                 cell_vectors,
                 color_dict,
                 marker_dict,
                 filename):
    """
    Function to plot elements change of a process dict.
    """
    title = process_dict["reaction"]
    coords = process_dict["coordinates_group"][0]
    elements_before = process_dict["elements_before"]
    elements_after = process_dict["elements_after"]
    basis = process_dict["basis_sites"][0]
    move_vector = basis_coords[basis]

    # Before.
    coords_before = CoordsGroup(coords, elements_before).move(move_vector)
    fig = plot_elements(cell_vectors,
                        coords_before,
                        color_dict,
                        marker_dict,
                        title=title+"(before)")
    fig.savefig(filename + "_0.png")
    plt.close(fig)

    # After.
    coords_after = CoordsGroup(coords, elements_after).move(move_vector)
    fig = plot_elements(cell_vectors,
                        coords_after,
                        color_dict,
                        marker_dict,
                        title=title+"(after)")
    fig.savefig(filename + "_1.png")
    plt.close(fig)


if __name__ == "__main__":
    cell_vectors = [[1.0, 0.0, 0.0],
                    [0.5, 0.87, 0.0],
                    [0.0, 0.0, 1.0]]

    marker_dict = dict(V='o',
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

    basis_coords = [[0.0, 0.0, 0.0],
                    [0.0, 0.5, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.5, 0.5, 0.0],
                    [1.0/3.0, 1.0/3.0, 0.0],
                    [2.0/3.0, 2.0/3.0, 0.0]]

    globs, locs = {}, {}
    execfile("kmc_processes.py", globs, locs)

    path = "./processes/"
    if not os.path.exists(path):
        os.mkdir(path)

    for idx, process in enumerate(locs["processes"]):
        filename = "{}process_{}".format(path, idx)
        logging.info("plotting process {}".format(process["reaction"]))
        plot_process(process, basis_coords, cell_vectors,
                     color_dict, marker_dict, filename)
        logging.info("Ok.")

