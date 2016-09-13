"""
Script to run kynetix kMC simulation.
"""

import commands
import logging
import time

from kynetix import mpi_master
from kynetix.model import KineticModel
from kynetix.utilities.format_utilities import convert_time

if "__main__" == __name__:
    # Remove old log file.
    commands.getstatusoutput("rm -rm out.log auto_*")

    # Set logger.
    logger = logging.getLogger("model.KMCModelRun")

    # Construct KMC model.
    model = KineticModel(setup_file="pt-111.mkm")
    parser = model.parser()
    parser.parse_data(filename="rel_energy.py", relative=True)

    start = time.time()
    try:
        model.run_kmc()
    except Exception as e:
        # Log exception info.
        if mpi_master:
            msg = "{} exception is catched.".format(type(e).__name__)
            logger.exception(msg)
        raise e

    end = time.time()
    t = end - start
    h, m, s = convert_time(t)

    if mpi_master:
        logger.info("Time used: {:d} h {:d} min {:f} sec".format(h, m, s))

