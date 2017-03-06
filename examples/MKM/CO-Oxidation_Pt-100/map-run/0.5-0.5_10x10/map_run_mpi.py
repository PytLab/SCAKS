#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands
import logging
import os
import time

import numpy as np
from mpi4py import MPI

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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    pO2 = np.linspace(1e-5, 0.5, 10)
else:
    pO2 = None
pO2 = comm.scatter(pO2, root=0)

pCOs = np.linspace(1e-5, 0.5, 10)
#pO2s = np.arange(0.01, 1.0, 0.02)

if "__main__" == __name__:
    # Clean up current dir.
    commands.getstatusoutput("rm -rf out.log auto_*")

    if rank == 0:
        start = time.time()

    tofs_2d = []
    try:
        # Construct setup dict.
        setup_dict['species_definitions']['O2_g']['pressure'] = pO2

        tofs_1d = []
        for j, pCO in enumerate(pCOs):
            print "INFO: runing pO2: {} pCO: {}".format(pO2, pCO)
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

        comm.gather(tofs_1d, root=0)

        if rank == 0:
            end = time.time()
            print " "

    finally:
        if rank == 0:
            # Write tofs to file.
            p_str = "pCO = {}\n\npO2 = {}\n\n".format(pCOs.tolist(), pO2s.tolist())
            adsorbates_str = "adsorbates = {}\n\n".format(model.adsorbate_names)
            essential_str = p_str + adsorbates_str

            tof_str = "tofs = {}\n\n".format(tofs_2d)
            with open("auto_tofs.py", "w") as f:
                f.write(essential_str + tof_str)

            delta_time = end - start
            h, m, s = convert_time(delta_time)
            print "Time used: {:d} h {:d} min {:f} sec ({:.2f}s)".format(h, m, s, delta_time)

