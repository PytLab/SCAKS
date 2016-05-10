'''
   Module containing error handling.
'''

from exceptions import Exception


class Error(Exception):
    """ Class for describing an error. """

    def __init__(self, msg=''):
        '''
        Constructor for the error class.

        :param msg: A message to the user describing what went wrong. If none
                    is given the message is set to an empty string.
        :type msg: string
        '''
        self.__msg = msg

    def __str__(self):
        """
        Get a string representation of the error.

        :returns: The error message string.
        """
        return self.__msg


# kinetic model errors
class ParameterError(Error):
    ''' Class for model paramter error. '''
    pass


class ToolsImportError(Error):
    ''' Class for errors in instantiation of model tools. '''
    pass


class SetupError(Error):
    ''' Class for type error of inputs in setup file. '''
    pass


# parser errors
class ProcessParsingError(Error):
    ''' Class for errors in parsing kMC process/event. '''
    pass


class SpeciesError(Error):
    ''' Class for species error in reaction equation. '''
    pass


class ElementSearchingError(Error):
    ''' Class for errors of element info searching in database. '''
    pass


class ReactionEquationError(Error):
    ''' Class for errors in reaction equations.'''
    pass


# solvers errors
class GridTypeError(Error):
    ''' Class for errors in grid type of lattice. '''
    pass


class FilesError(Error):
    ''' Class for errors in files. '''
    pass
