"""
Script for running Micro-kinetic Model simulation.
"""

import commands
import logging
import sys
import time

from kynetix import mpi_master
from kynetix.model import KineticModel
from kynetix.utilities.format_utilities import convert_time

# Custom parameters.
UseRelativeEnergy = True    # Use only relative energies.
OdeEnd = 10000              # ODE integration time limit.
CalcXRC = False             # Calculate Degree of Rate Control(XRC) or not.
ProductionName = "CH3OH_g"  # Production name of your model.
OdeOnly = False             # Do ODE integration only.

if "__main__" == __name__:
    # Clean up current dir.
    commands.getstatusoutput("rm -rf out.log data.pkl auto_*")

    # Set script logger.
    logger = logging.getLogger("model.MkmRunScript")

    # Get setup file.
    status, output= commands.getstatusoutput("ls *.mkm | tail -1")
    if status:
        if mpi_master:
            logger.error(output)
            logger.info("Exiting...")
        sys.exit(1)

    start = time.time()
    try:
        # Build micor-kinetic model.
        model = KineticModel(setup_file=output)

        # Read data.
        parser = model.parser()
        solver = model.solver()
        parser.parse_data(relative=UseRelativeEnergy)
        solver.get_data()

        # Initial coverages guess.
        trajectory = solver.solve_ode(time_end=OdeEnd)
        init_guess = trajectory[-1]

        # Run.
        model.run_mkm(init_cvgs=init_guess,
                      solve_ode=OdeOnly,
                      coarse_guess=False,
                      relative=True,
                      XRC=CalcXRC,
                      product_name=ProductionName)
    except Exception as e:
        if mpi_master:
            msg = "{} exception is catched.".format(type(e).__name__)
            logger.exception(msg)
        raise e

    # Time used.
    end = time.time()
    t = end - start
    h, m, s = convert_time(t)

    if mpi_master:
        logger.info("Time used: {:d} h {:d} min {:f} sec".format(h, m, s))

