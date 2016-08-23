import time
try:
    import cPickle as pickle
except ImportError:
    import picke

try:
    from mpi4py import MPI
    mpi_installed = True
    mpi_comm = MPI.COMM_WORLD
    mpi_rank = mpi_comm.Get_rank()
    mpi_size = mpi_comm.Get_size()
except ImportError:
    mpi_installed = False
    mpi_rank = 0
    mpi_size = 1

from kynetix.functions import *
from kynetix.errors.error import *


# Condition for info output or not.
mpi_master = (mpi_rank == 0)

__version__ = '1.0.0'

file_header = (
    '# This file was automatically generated by Kynetix' +
    ' (https://github.com/PytLab/Kynetix).\n' +
    '# Version %s\n# Date: %s \n#\n' +
    '# Do not make changes to this file ' +
    'unless you know what you are doing\n\n') % (__version__, time.asctime())

#-------------------------------------------------------
# Some base classes for kinetic model are defined below |
#-------------------------------------------------------


class ModelShell(object):
    """
    A non-functional parent class to be inherited by
    other tools class of kinetic model.
    """

    def __init__(self, owner):
        self._owner = owner
        self._archived_data_dict = {}

    def update_defaults(self, defaults):
        """
        Update values in defaults dict,
        if there are custom parameters in setup file.

        Parameters:
        -----------
        default: default attributes dict, dict.
        """

        for parameter_name in defaults:
            attribute_name = mangled_name(self._owner, parameter_name)
            if hasattr(self._owner, attribute_name):
                defaults[parameter_name] = getattr(self._owner, attribute_name)

        return defaults

    def archive_data(self, data_name, data):
        """
        Update data dict and dump it to data file.

        Parameters:
        -----------
        data_name: key in data dict, str.

        data: value in data dict, any python data type.
        """
        # Update data dict.
        if data_name in self._archived_variables:
            self._archived_data_dict[data_name] = data
            # Dump data dict to data file
            if self._archived_data_dict:
                with open(self._owner.data_file(), 'wb') as f:
                    cPickle.dump(self._archived_data_dict, f)

    @staticmethod
    def write2file(filename, line):
        f = open(filename, 'a')
        f.write(line)
        f.close()

