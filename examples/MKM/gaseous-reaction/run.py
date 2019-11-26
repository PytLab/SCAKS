"""
Script for running Micro-kinetic Model simulation.
"""
import pickle as pkl

from scaks.models.micro_kinetic_model import MicroKineticModel
from scaks.solvers.steady_state_solver import SteadyStateSolver
from scaks.correctors.thermodynamic_corrector import ThermodynamicCorrector

model_dict = {
    'rxn_expressions': [
        'H2O_g + *_o + *_b <-> H-OH_b + *_o -> H_o + OH_b',
        'O2_g + *_n -> O2_n',
        'O2_n + CH4_g + *_b <-> CH3-H_g + O2_n + *_b -> CH3_i + O_n + OH_b',
        'CH3_i + O_n <-> CH3-O_n -> OCH3_n',
        'OCH3_n + OH_b + H_o <-> OCH3_n + OH-H_b + *_o -> OCH3_i + H2O_g + *_n + *_b + *_o',
        'CH4_g + OH_b <-> CH3-H_g + OH_b -> CH3_i + H2O_g + *_b',
        'OCH3_i <-> OCH2-H_i -> HCHO_g + H_i', #1
        'OCH3_i + O2_g <-> OCH2-H_i + O2_g -> HCHO_g + OOH_i', #2
        'HCHO_g + CH3_i <-> CH3_i + H-CHO_g -> CHO_i + CH4_g', #3
        'CHO_i + O2_g <-> CO-H_i + O2_g -> CO_g + OOH_i', #4
        'CHO_i <-> CO-H_i -> CO_g + H_i', #5
        '2CH3_i -> C2H6_g', #6
        '2H_i -> H2_g', #7
        'H_i + CH4_g <-> CH3-H_g + H_i -> H2_g + CH3_i', #8
        'H_i + C2H6_g <-> C2H5-H_g + H_i -> H2_g + C2H5_i', #9
        'O2_g + H_i <-> OO-H_i -> OOH_i', #10
        'CH3_i + OOH_i -> CH3OOH_g', #11
        'CH3OOH_g <-> CH3O-OH_g -> CH3O_i + OH_i', #12
        'OH_i + H_i -> H2O_g', #13
        'C2H6_g + CH3_i <-> C2H5-H_g + CH3_i -> C2H5_i + CH4_g', #14
        'C2H5_i -> C2H4_g + H_i', #15
    ],

    'species_definitions': {
        'CH4_g': {'pressure': 0.28},
        'O2_g': {'pressure': 0.14},
        'CO_g': {'pressure': 0.043},
        'H2_g': {'pressure': 0.0078},
        'C2H4_g': {'pressure': 0.0054},
        'C2H6_g': {'pressure': 0.0054},
        'CO2_g': {'pressure': 0.0025},
        'H2O_g': {'pressure': 0.0},
        'HCHO_g': {'pressure': 0.0},
        'CH3OOH_g': {'pressure': 0.0},
        '*_b': {'site_name': 'B', 'type': 'site', 'total': 1.0},
        '*_o': {'site_name': 'O', 'type': 'site', 'total': 1.0},
        '*_n': {'site_name': 'N', 'type': 'site', 'total': 1.0},
    },

    'temperature': 690 + 273.15,
    'parser': "RelativeEnergyParser",

    'rate_algo': 'CT',
    'unitcell_area': 9.0e-20,
    'active_ratio': 1.0,

    'ode_output_interval': 100,
    'decimal_precision': 1000,
    'tolerance': 1e-15
}

# Build micor-kinetic model.
model = MicroKineticModel(setup_dict=model_dict)

# Read data.
model.parser.parse_data(filename='./rel_energy.py')

# Create solver
solver = SteadyStateSolver(model)
model.set_solver(solver)
model.solver.get_data()

# Create corrector for energy correction
corrector = ThermodynamicCorrector(model)
model.set_corrector(corrector)

#model.solver.get_data_symbols()

def ss_equation(t, y):
    return model.solver.steady_state_function(y)

if __name__ == '__main__':

    cvgs = [0.0]*13
    trajectory = model.solver.solve_ode(algo='lsoda', time_span=0.1, initial_cvgs=cvgs, time_end=100, traj_output=True)
    init_guess = trajectory[-1]

    # Run.
    model.run(init_cvgs=init_guess,
              solve_ode=False,
              coarse_guess=False,
              XRC=False,
              product_name='CO_g')

