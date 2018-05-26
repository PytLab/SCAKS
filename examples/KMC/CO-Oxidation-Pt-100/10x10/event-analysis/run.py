"""
Script to run kynetix kMC simulation.
"""

import logging
import time

from kynetix.compatutil import subprocess
from kynetix.mpicommons import mpi
from kynetix.models.kmc_model import KMCModel
from kynetix.utilities.format_utilities import convert_time

if "__main__" == __name__:
    # Remove old log file.
    subprocess.getstatusoutput("rm -rm out.log auto_*")

    # Set logger.
    logger = logging.getLogger("model.KMCModelRun")

    # Construct KMC model.
    model = KMCModel(setup_file="pt-100.mkm")
    parser = model.parser
    parser.parse_data(energy_file="rel_energy.py")

    start = time.time()
    try:
        model.run()
    except Exception as e:
        # Log exception info.
        if mpi.is_master:
            msg = "{} exception is catched.".format(type(e).__name__)
            logger.exception(msg)
        raise e

    end = time.time()
    t = end - start
    h, m, s = convert_time(t)

    if mpi.is_master:
        logger.info("Time used: {:d} h {:d} min {:f} sec".format(h, m, s))

