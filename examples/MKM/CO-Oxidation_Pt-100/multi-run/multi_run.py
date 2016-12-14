#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os

import numpy as np

from kynetix.models.micro_kinetic_model import MicroKineticModel

setup_dict = dict(
    rxn_expressions = [
        'CO_g + *_s -> CO_s',
        'O2_g + 2*_s -> O2_2s',
        'O2_2s + CO_s <-> OCO-O_2s + *_s -> O_s + CO2_g + 2*_s',
        'O2_g + 2*_s -> 2O_s',
        'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
    ],

    species_definitions = {
        'CO_g': {'pressure': 0.10},
        'O2_g': {'pressure': 0.2},
        'CO2_g': {'pressure': 0.0288},
        's': {'site_name': 'top', 'type': 'site', 'total': 1.0},
    },

    temperature = 500,

    unitcell_area = 9.0e-20,
    active_ratio = 4./9.,

    parser = "RelativeEnergyParser",
    solver = "SteadyStateSolver",
    corrector = "ThermodynamicCorrector",
    plotter = "EnergyProfilePlotter",

    rate_algo = "CT",
    rootfinding = "MDNewton",
    tolerance = 1e-15,
    max_rootfinding_iterations = 100,
)

pCOs = np.arange(0.01, 0.2, 0.001)

if "__main__" == __name__:
    logger = logging.getLogger("model.MultiRun")

    ss_cvgs = []
    tofs = []
    for i, pCO in enumerate(pCOs):
        # Construct setup dict.
        setup_dict['species_definitions']['CO_g']['pressure'] = pCO

        # Construct model.
        model = MicroKineticModel(setup_dict=setup_dict, verbosity=logging.INFO)

        # Read data.
        model.parser.parse_data(relative=True)
        model.solver.get_data()

        # Initial coverage guess.
        if i == 0:
            trajectory = model.solver.solve_ode(time_span=0.0001,
                                                time_end=10,
                                                traj_output=False)
            init_guess = trajectory[-1]
        else:
            init_guess = ss_cvgs[-1]

        # Run.
        logger.info("Running pressure CO_g: {}".format(pCO))
        model.run(init_cvgs=init_guess,
                  solve_ode=False,
                  relative=True,
                  XRC=False,
                  product_name="CO2_g")
        model.clear_handlers()

        # Collect TOF.
        tof_idx = model.gas_names.index("CO2_g")
        tofs.append(float(model.TOFs[tof_idx]))

        # Collect steady state coverages.
        ss_cvgs.append(model.steady_state_coverages)

    # Write tofs to file.
    tof_str = "tofs = {}\n".format(tofs)
    p_str = "pCO = {}\n".format(pCOs.tolist())
    with open("auto_tofs.py", "w") as f:
        f.write(tof_str + p_str)

