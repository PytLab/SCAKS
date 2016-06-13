import logging

from kynetix.errors.error import *
from kynetix.plotters.en_profile import *
from kynetix.plotters.plotter_base import *


class EnergyProfilePlotter(PlotterBase):
    def __init__(self, owner):
        super(EnergyProfilePlotter, self).__init__(owner)

        # Set logger.
        self.__logger = logging.getLogger("model.plotters.KMCPlotter")

    def plot_single(self, rxn_expression, **kwargs):
        """
        Function to plot energy profile for an elementary reaction.

        Parameters:
        -----------
        rxn_expression: An elementary reaction epxression string, str.
                        Must be one of owner's rxn_expressions.

        kwargs: See doc string of kynetix.plotters.en_profile.plot_single_energy_diagram.

        Returns:
        --------
        Matplotlib.figure.Figure object.
        """
        # {{{
        rxn_expressions = self._owner.rxn_expressions()

        # Check.
        if rxn_expression not in self._owner.rxn_expressions():
            msg = "'{}' not in model's reaction expressions.".format(rxn_expression)
            raise ParameterError(msg)

        # Get energy tuple.
        idx = rxn_expressions.index(rxn_expression)

        if not self._owner.has_relative_energy():
            msg_template = "Model '{}' has no relative energy, please try to parse data."
            msg = msg.format(self._owner)
            raise AttributeError(msg)

        # Convert relative energies to energy tuple.
        relative_energies = self._owner.relative_energies()
        Gaf = relative_energies['Gaf'][idx]
        dG = relative_energies['dG'][idx]
        energy_tuple = (0.0, Gaf, dG) if Gaf else (0.0, dG)

        # Plot.
        fig, _, _ = plot_single_energy_diagram(energy_tuple, rxn_expression, **kwargs)

        return fig
        # }}}

    def plot_all(self, **kwargs):
        """
        Function to plot a merged energy profile with all elementary reactions.

        Parameters:
        -----------
        kwargs: See doc string of kynetix.plotters.en_profile.plot_single_energy_diagram.

        Returns:
        --------
        Matplotlib.figure.Figure object.
        """
        # Get reaction expressions.
        rxn_expressions = self._owner.rxn_expressions()

        # Get relative energies.
        if not self._owner.has_relative_energy():
            msg_template = "Model '{}' has no relative energy, please try to parse data."
            msg = msg.format(self._owner)
            raise AttributeError(msg)

        relative_energies = self._owner.relative_energies()
        energy_tuples = []
        # Loop over all relative energies to collect energy tuples.
        for Gaf, dG in zip(relative_energies['Gaf'], relative_energies['dG']):
            energy_tuple = (0.0, Gaf, dG) if Gaf else (0.0, dG)
            energy_tuples.append(energy_tuple)

        # Plot.
        fig, _, _ = plot_multi_energy_diagram(rxn_expressions, energy_tuples, **kwargs)

        return fig

