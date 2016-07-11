from kynetix import ModelShell


class PlotterBase(ModelShell):
    def __init__(self, owner):
        """
        A class acts as a base class to be inherited by other
        plotter classes, it is not functional on its own.
        """
        ModelShell.__init__(self, owner)

