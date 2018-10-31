import logging
import os

from ..parsers.rxn_parser import RxnEquation
from ..errors.error import *
from .plotter_base import *
from catplot.ep_components.ep_canvas import EPCanvas
from catplot.ep_components.ep_lines import ElementaryLine


class EnergyProfilePlotter(PlotterBase):
    def __init__(self, owner):
        ''' Class for energy profile plotting

        :param owner: The kinetic model that own this plotter
        :type owner: KineticModel
        '''
        super(EnergyProfilePlotter, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger("model.plotters.EnergyProfilePlotter")

    def plot_all(self, **kwargs):
        ''' Plot energy profile for all elementary reactions.

        :param energies: energies for states of a elementary reaction
        :type energies: tuple or list

        :param n: the point number in each state, default is 100.
        :type n: int

        :param hline_leng: the length of the horizontal line for the IS and FS
        :type hline_length: float

        :param peak_width: the width of the peak in energy profile, default is 1.0.
        :type peak_width: float

        :param interp_method: the type of interpolation algorithm, possible value: "spline", "quadratic", default is "spline".
        :type interp_method: str

        :param rxn_equation: elementary reaction equation, default is None.
        :type rxn_equation: str

        :param line_width: line width, default is 3.
        :type line_width: float

        :param color: color code of the line, default is #000000 (black).
        :type color: str

        :param shadow_color: color code of the shadow lines, default is #595959.
        :type shadow_color: str

        :param shadow_depth: shadow depth of the line, default is 0, no shadow.
        :type shadow_depth: int
        '''
        if not os.path.exists('energy_profile'):
            os.mkdir('energy_profile')

        for idx, (rxn, Ga, dG) in enumerate(zip(self._owner.rxn_expressions,
                                                self._owner.relative_energies['Gaf'],
                                                self._owner.relative_energies['dG'])):
            energies = [0.0, dG] if abs(Ga) < 1e-5 else [0.0, Ga, dG]
            canvas = EPCanvas()
            canvas.add_line(ElementaryLine(energies, rxn_equation=rxn, **kwargs))
            canvas.draw()
            image_name = 'energy_profile/{}_{}.png'.format(idx, rxn)
            canvas.figure.savefig(image_name)
            self.__logger.info('{} is saved'.format(image_name))

