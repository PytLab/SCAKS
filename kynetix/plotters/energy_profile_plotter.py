from kynetix.plotters.en_profile import *
from kynetix.plotters.plotter_base import *


class EnergyProfilePlotter(PlotterBase):
    def __init__(self, owner):
        super(EnergyProfilePlotter, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger("model.plotters.KMCPlotter")

