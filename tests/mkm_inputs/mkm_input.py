#!/usr/bin/env python
# -*- coding: utf-8 -*-

rxn_expressions = [
    'CO_g + *_s -> CO_s',
    'O2_g + 2*_s -> 2O_s',
    'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
]

species_definitions = {
    'CO_g': {'pressure': 1.0},
    'O2_g': {'pressure': 1./3.},
    'CO2_g': {'pressure': 0.00},
    '*_s': {'site_name': '111', 'type': 'site', 'total': 1.0},
}

temperature = 450.0
parser = "AbsoluteEnergyParser"
solver = "SteadyStateSolver"
corrector = "ThermodynamicCorrector"
plotter = "EnergyProfilePlotter"
rootfinding = 'ConstrainedNewton'
decimal_precision = 100
tolerance = 1e-20
max_rootfinding_iterations = 100

