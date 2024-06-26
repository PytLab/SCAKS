""" Module providing interface for profiling.
"""

import cProfile
import pstats
import os

def do_cprofile(filename):
    """ Decorator for function profiling.
    """
    def wrapper(func):
        def profiled_func(*args, **kwargs):

            # Flag for do profiling or not.
            DO_PROF = os.getenv("PROFILING")

            if DO_PROF:
                profile = cProfile.Profile()
                profile.enable()
                result = func(*args, **kwargs)
                profile.disable()
                # Sort stat by internal time.
                sortby = "tottime"
                ps = pstats.Stats(profile).sort_stats(sortby)
                ps.dump_stats(filename)
            else:
                result = func(*args, **kwargs)
            return result
        return profiled_func
    return wrapper


def print_pstats(filename, sortby='cumtime'):
    """ Function to print profiling stats results.
    """
    p = pstats.Stats(filename)
    p.sort_stats(sortby).print_stats(10, 1.0, '.*')

