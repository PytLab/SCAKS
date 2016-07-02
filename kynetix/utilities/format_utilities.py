"""
Module for holding common checking utility functions.
"""

import logging

from kynetix.errors.error import *


def get_list_string(var_name, list_obj, ncols=5):
    '''
    Function to get string format of a list object.

    Parameters:
    -----------
    var_name: name of variable, str.

    list_obj: object of variable.

    ncols: number of columns. 

    Example:
    --------
    >>> a = [(1,2,3), (2,3,4)]
    >>> get_list_string('a', a)
    >>> 'a = [\n    (1, 2, 3),\n    (2, 3, 4),\n]\n\n'

    Called:
    -------
    kmc_plugins.CoveragesAnalysis()

    '''
    begin = var_name + ' = ['
    indent = ' '*4
    data = ''
    for idx, elem in enumerate(list_obj):
        # if item is iterable, one a line
        if hasattr(elem, '__iter__'):
            data += ('\n' + indent + str(elem) + ',')
            continue
        # 5 items a line by default
        if idx % ncols == 0:
            data += ('\n' + indent)

        # add single quotes for string
        if isinstance(elem, str):
            data += ("'" + elem + "', ")
        else:
            data += (str(elem) + ', ')
    end = '\n]\n\n'

    content = begin + data + end

    return content

def get_dict_string(var_name, dict_obj):
    """
    Function to get string format of a list object.

    Parameters:
    -----------
    var_name: name of variable, str.

    dict_obj: the dict object whose string would be returned.
    """
    begin = var_name + " = {\n"
    data = ""
    for key in sorted(dict_obj):
        value = dict_obj[key]
        data += "    '{}': {},\n".format(key, value)
    end = "}\n\n"

    content = begin + data + end

    return content


def convert_time(sec):
    """
    Convert format of time from seconds to *h *min *sec.
    """
    hours = int(sec/(3600.0))
    minutes = int((sec - hours*3600)/60.0)
    seconds = sec - hours*3600 - minutes*60

    # int, int, float
    return hours, minutes, seconds

