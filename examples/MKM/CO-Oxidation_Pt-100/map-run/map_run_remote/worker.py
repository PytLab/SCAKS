#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands
import logging
import os
import time
import Queue
from multiprocessing import Pool
from multiprocessing.managers import BaseManager

import numpy as np

from kynetix.models.micro_kinetic_model import MicroKineticModel

ADDR = '127.0.0.1'
PORT = 5000
AUTHKEY = 'pytlab'

setup_dict = dict(
    # {{{
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
    # }}}
)

#pCOs = np.linspace(1e-5, 0.5, 10)
#pO2s = np.linspace(1e-5, 0.5, 10)

def get_manager():
    class WorkManager(BaseManager):
        pass

    WorkManager.register('get_jobid_queue')
    WorkManager.register('get_tofs_list')
    WorkManager.register('get_pCOs')
    WorkManager.register('get_pO2s')

    manager = WorkManager(address=(ADDR, PORT), authkey=AUTHKEY)

    return manager

def task(pO2):
    # {{{
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
        model.run(init_cvgs=init_guess, product_name="CO2_g")
        model.clear_handlers()

        # Collect TOF.
        tof_idx = model.gas_names.index("CO2_g")
        tofs_1d.append(float(model.TOFs[tof_idx]))

    return tofs_1d
    # }}}

if "__main__" == __name__:

    manager = get_manager()
    print "work manager connect to {}:{}...".format(ADDR, PORT)
    manager.connect()

    shared_tofs_list = manager.get_tofs_list()
    shared_jobid_queue = manager.get_jobid_queue()
    pCOs = manager.get_pCOs()
    shared_pO2s = manager.get_pO2s()

    pool = Pool()

    while 1:
        try:
            indices = shared_jobid_queue.get_nowait()
            pO2s = [shared_pO2s[i] for i in indices]
            print "Run {}".format(str(pO2s))
#            import ipdb
#            ipdb.set_trace()
            tofs_2d = pool.map(task, pO2s)

            # Update shared tofs list.
            for idx, tofs_1d in zip(indices, tofs_2d):
                shared_tofs_list[idx] = tofs_1d
        except Queue.Empty:
            break
    manager.shutdown()

