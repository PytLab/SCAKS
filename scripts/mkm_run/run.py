import commands
import logging
import sys
import time

from kynetix import mpi_master
from kynetix.model import KineticModel
from kynetix.utilities.format_utilities import convert_time

# Custom parameters.
UseRelativeEnergy = True
OdeEnd = 10000

if "__main__" == __name__:
    # Clean up current dir.
    commands.getstatusoutput("rm -rf out.log data.pkl auto_*")

    # Set script logger.
    logger = logging.getLogger("model.MkmRun")

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
        model = KineticModel(setup_file=setup_file)

        # Read data.
        parser = model.parser()
        solver = model.solver()
        parser.parse_data(relative=UseRelativeEnergy)
        solver.get_data()

        # Initial coverages guess.
        ss = solver.solve_ode(time_end=OdeEnd)
        c = ss[-1]

        # Run.
        model.run_mkm(init_cvgs=cvgs, coarse_guess=False, relative=True)
    except Exception as e:
        if mpi_master:
            msg = "{} exception is catched.".format(type(e).__name__)
            logger.exception(msg)
        raise e

    end = time.time()
    t = end - start
    h, m, s = convert_time(t)

    if mpi_master:
        logger.info("Time used: {:d} h {:d} min {:f} sec".format(h, m, s))

