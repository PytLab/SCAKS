import logging
import time

from kynetix.model import KineticModel
from kynetix.utilities.format_utilities import convert_time

if "__main__" == __name__:

    # Set logger.
    logger = logging.getLogger("ModelRun")

    # Construct KMC model.
    model = KineticModel(setup_file="pt-100.mkm")
    parser = model.parser()
    parser.parse_data(filename="rel_energy.py", relative=True)

    start = time.time()
    try:
        model.run_kmc()
    except Exception as e:
        # Log exception info.
        msg = "{} exception is catched.".format(type(e).__name__)
        logger.exception(msg)
        raise e

    end = time.time()
    t = end - start
    h, m, s = convert_time(t)

    logger.info("Time used: {:d} h {:d} min {:f} sec".format(h, m, s))

