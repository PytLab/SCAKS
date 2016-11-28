"""
Module providing interface for profiling.
"""

import cProfile
import pstats
import os

def do_cprofile(filename):
    """
    Decorator for function profiling.
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
                # Sort stat by cumulative time.
                sortby = "cumulative"
                ps = pstats.Stats(profile).sort_stats(sortby)
                ps.dump_stats(filename)
                ps.print_stats()
            else:
                result = func(*args, **kwargs)

            return result
        return profiled_func
    return wrapper

