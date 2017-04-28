#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import time

import numpy as np

from kynetix.compatutil import subprocess
from kynetix.models.micro_kinetic_model import MicroKineticModel
from kynetix.utilities.format_utilities import convert_time

setup_dict = dict(
    rxn_expressions = [
        'CO_g + *_s -> CO_s',
#        'O2_g + 2*_s -> O2_2s',
#        'O2_2s + CO_s <-> OCO-O_2s + *_s -> O_s + CO2_g + 2*_s',
        'O2_g + 2*_s -> 2O_s',
        'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
    ],

    species_definitions = {
        'CO_g': {'pressure': 0.10},
        'O2_g': {'pressure': 0.2},
        'CO2_g': {'pressure': 0.01},
        '*_s': {'site_name': 'top', 'type': 'site', 'total': 1.0},
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

pCOs = np.linspace(1e-5, 0.5, 10)
pO2s = np.linspace(1e-5, 0.5, 10)
#pO2s = np.arange(0.01, 1.0, 0.02)

if "__main__" == __name__:
    # Clean up current dir.
    subprocess.getstatusoutput("rm -rf out.log auto_*")

    start = time.time()

    cvgs_CO_2d = []
    cvgs_O_2d = []
    tofs_2d = []
    try:
        for i, pO2 in enumerate(pO2s):
            # Construct setup dict.
            setup_dict['species_definitions']['O2_g']['pressure'] = pO2

            cvgs_CO_1d = []
            cvgs_O_1d = []
            tofs_1d = []
            for j, pCO in enumerate(pCOs):
                print("INFO: runing pO2: {} pCO: {}".format(pO2, pCO))
                setup_dict['species_definitions']['CO_g']['pressure'] = pCO

                # Construct model.
                model = MicroKineticModel(setup_dict=setup_dict,
                                          console_handler_level=logging.WARNING)

                # Read data.
                model.parser.parse_data()
                model.solver.get_data()

                # Initial coverage guess.
                trajectory = model.solver.solve_ode(time_span=0.0001,
                                                    time_end=10,
                                                    traj_output=False)
                init_guess = trajectory[-1]

                # Run.
                model.run(init_cvgs=init_guess,
                          product_name="CO2_g")
                model.clear_handlers()

                # Collect TOF.
                tof_idx = model.gas_names.index("CO2_g")
                tofs_1d.append(float(model.TOFs[tof_idx]))

                # Collect CO_s coverage.
                cvg_CO = float(model.steady_state_coverages[0])
                cvgs_CO_1d.append(cvg_CO)

                # Collect O_s coverage.
                cvg_O = float(model.steady_state_coverages[1])
                cvgs_O_1d.append(cvg_O)

            tofs_2d.append(tofs_1d)
            cvgs_CO_2d.append(cvgs_CO_1d)
            cvgs_O_2d.append(cvgs_O_1d)

            end = time.time()
            print(" ")

    finally:
        # Write tofs to file.
        p_str = "pCO = {}\n\npO2 = {}\n\n".format(pCOs.tolist(), pO2s.tolist())
        adsorbates_str = "adsorbates = {}\n\n".format(model.adsorbate_names)
        essential_str = p_str + adsorbates_str

        tof_str = "tofs = {}\n\n".format(tofs_2d)
        with open("auto_tofs.py", "w") as f:
            f.write(essential_str + tof_str)

        cvgs_O_str = "cvgs_O = {}\n\n".format(cvgs_O_2d)
        cvgs_CO_str = "cvgs_CO = {}\n\n".format(cvgs_CO_2d)
        with open("auto_cvgs.py", "w") as f:
            f.write(essential_str + cvgs_O_str + cvgs_CO_str)

        delta_time = end - start
        h, m, s = convert_time(delta_time)
        print("Time used: {:d} h {:d} min {:f} sec ({:.2f}s)".format(h, m, s, delta_time))

