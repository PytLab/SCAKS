#!/usr/bin/env python
# -*- coding: utf-8 -*-

import commands
import logging
import os
import time
from multiprocessing.managers import BaseManager, ListProxy
from multiprocessing import Queue

import numpy as np

from kynetix.models.micro_kinetic_model import MicroKineticModel


ADDR = ''
PORT = 5000
AUTHKEY = 'pytlab'
N = 10
NNODE = 2
pCOs = np.linspace(1e-5, 0.5, N).tolist()
pO2s = np.linspace(1e-5, 0.5, N).tolist()
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


def get_manager():
    class JobManager(BaseManager):
        pass

    jobid_queue = Queue()
    JobManager.register('get_jobid_queue', callable=lambda: jobid_queue)

    tofs = [None]*N
    JobManager.register('get_tofs_list', callable=lambda: tofs, proxytype=ListProxy)
    JobManager.register('get_pCOs', callable=lambda: pCOs, proxytype=ListProxy)
    JobManager.register('get_pO2s', callable=lambda: pCOs, proxytype=ListProxy)

    manager = JobManager(address=(ADDR, PORT), authkey=AUTHKEY)

    return manager

def fill_jobid_queue(manager, nclient):
    indices = range(N)
    interval = N/nclient
    jobid_queue = manager.get_jobid_queue()
    for i in range(nclient-1):
        jobid_queue.put(indices[i: i+interval])
    jobid_queue.put(indices[i+interval:])

def run_server():
    manager = get_manager()
    print "Start manager at {}:{}...".format(ADDR, PORT)
    manager.start()
    shared_job_queue = fill_jobid_queue(manager, NNODE)
    shared_tofs_list = manager.get_tofs_list()

    while None in shared_tofs_list:
        pass

    manager.shutdown()
    print "Manager shutdown"


if "__main__" == __name__:
    # Clean up current dir.
    commands.getstatusoutput("rm -rf out.log auto_*")

    start = time.time()
    run_server()
    end = time.time()
    delta_time = end - start

    model = MicroKineticModel(setup_dict=setup_dict,
                              console_handler_level=logging.WARNING)
    # Write tofs to file.
    p_str = "pCO = {}\n\npO2 = {}\n\n".format(pCOs.tolist(), pO2s.tolist())
    adsorbates_str = "adsorbates = {}\n\n".format(model.adsorbate_names)
    essential_str = p_str + adsorbates_str

    tof_str = "tofs = {}\n\n".format(tofs_2d)
    with open("auto_tofs.py", "w") as f:
        f.write(essential_str + tof_str)

    delta_time = end - start
    h, m, s = convert_time(delta_time)
    print "Time used: {:d} h {:d} min {:f} sec ({:.2f} s)".format(h, m, s, delta_time)

