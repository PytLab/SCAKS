import os, sys 
import inspect
from functions import *


__all__ = ['loggers', 'parsers', 'solvers']

__version__ = '0.0.1'


class ModelShell(object):
    """
    A non-functional parent class to be inherited by
    other tools class of kinetic model.
    """
    def __init__(self, owner):
        self._owner = owner
        self.set_logger()

    def set_logger(self):
        #create a logger
        #import logger and get an instance of Logger class
        #The 'BLACK MAGIC' is hacked from CatMap (catmap/model.py)
        #https://docs.python.org/2/library/functions.html#__import__
        basepath = os.path.dirname(
            inspect.getfile(inspect.currentframe()))
        if basepath not in sys.path:
            sys.path.append(basepath)
        #from loggers import logger
        _module = __import__('loggers.logger', globals(), locals())
        logger_instance = getattr(_module, 'Logger')(owner=self)
        setattr(self, 'logger', logger_instance)

    def split_species(self, species_str):
        "Split species_str to number(int) and species_name(str)"
        if not '*' in species_str:  # for species adsorbated on surface
            m = self._owner.regex_dict['species'][0].search(species_str)
            if not m.group(1):
                stoichiometry = 1
            else:
                stoichiometry = int(m.group(1))
            species_name = m.group(2) + '_' + m.group(3)
            return stoichiometry, species_name
        else:  # for site
            m = self._owner.regex_dict['empty_site'][0].search(species_str)
            if not m.group(1):
                stoichiometry = 1
            else:
                stoichiometry = int(m.group(1))
            site_name = '*_' + m.group(2)
            return stoichiometry, site_name

    def update_defaults(self, defaults):
        """
        Update values in defaults dict,
        if there are custom parameters in setup file.
        """
        for parameter_name in defaults:
            if hasattr(self._owner, parameter_name):
                defaults[parameter_name] = \
                    getattr(self._owner, parameter_name)

        return defaults
