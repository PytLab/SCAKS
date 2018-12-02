"""
Script for running Micro-kinetic Model simulation.
"""

import logging
import time

from scaks.mpicommons import mpi
from scaks.models.micro_kinetic_model import MicroKineticModel
from scaks.utilities.format_utilities import convert_time

# Set script logger.
logger = logging.getLogger("model.MkmRunScript")

# Build micor-kinetic model.
model = MicroKineticModel(setup_file='model.py')

# Read data.
model.parser.parse_data()

# Pass data to solver
model.solver.get_data()

if __name__ == '__main__':
    # Initial coverages guess using ODE integration
    trajectory = model.solver.solve_ode(time_span=0.01, time_end=1, traj_output=True)
    init_guess = trajectory[-1]

    start = time.time()

    # Run.
    model.run(init_cvgs=init_guess,
              solve_ode=False,
              coarse_guess=False,
              XRC=False,
              product_name='CO2_g')

    # Time used.
    end = time.time()
    t = end - start
    h, m, s = convert_time(t)

    logger.info("Time used: {:d} h {:d} min {:f} sec".format(h, m, s))

