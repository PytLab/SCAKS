from pynetics import ModelShell


class CorrectorBase(ModelShell):
    def __init__(self, owner):
        """
        A class acts as a base class to be inherited by other
        corrector classes, it is not functional on its own.
        """
        ModelShell.__init__(self, owner)
