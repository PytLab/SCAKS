rxn_expressions = [
    'CO_g + *_s -> CO_s',
    'O2_g + 2*_s <-> O-O_2s -> 2O_s',
    'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s',
]

# Gas pressure.
species_definitions = {}
species_definitions['CO_g'] = {'pressure': 1.32*10**-7}
species_definitions['O2_g'] = {'pressure': 5.26*10**-7}
species_definitions['CO2_g'] = {'pressure': 1.32*10**-7}

# Site info.
species_definitions['*_s'] = {'site_name': 'top', 'type': 'site', 'total': 1.0}

# Temperature.
temperature = 500  # K

unitcell_area = 9.0e-20
active_ratio = 4./9.

parser = "RelativeEnergyParser"
solver = "SteadyStateSolver"
corrector = "ThermodynamicCorrector"
plotter = "EnergyProfilePlotter"

rate_algo = "CT"
rootfinding = "MDNewton"
tolerance = 1e-50
max_rootfinding_iterations = 100

