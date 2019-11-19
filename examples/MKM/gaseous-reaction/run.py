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
        'O2_n + CH4_g + *_b <-> CH3-H_g + O2_n -> CH3_i + O_n + OH_b',
        'CH3_i + O_n <-> CH3-O_n -> OCH3_n',
    ],

    'species_definitions': {
        'H2_g': {'pressure': 2.},
        'C8H6NO2_g': {'pressure': 0.005},
        'H2O_g': {'pressure': 0.01},
        'C8H6NH2_g': {'pressure': 0.005},
        'C8H6NHOH_g': {'pressure': 0.1},
        '*_a': {'site_name': 'TiO2', 'type': 'site', 'total': 1.0},
        '*_b': {'site_name': 'interface', 'type': 'site', 'total': 1.0},
        '*_c': {'site_name': 'Au', 'type': 'site', 'total': 1.0},
    },

    'temperature': 353.15,
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

model.solver.get_data_symbols()

def ss_equation(t, y):
    return model.solver.steady_state_function(y)

if __name__ == '__main__':

    from results import cvgs
#    with open('root.pkl', 'rb') as f:
#        data = pkl.load(f)
#    init_guess = cvgs = data['steady_state_coverages']
    trajectory = model.solver.solve_ode(algo='lsoda', time_span=0.1, initial_cvgs=cvgs, time_end=10**3, traj_output=True)
    init_guess = trajectory[-1]

    # Run.
    model.run(init_cvgs=init_guess,
              solve_ode=False,
              coarse_guess=False,
              XRC=False,
              product_name='C8H6NH2_g')
    #model.solver.get_GXRC(19, cvgs, model.transition_state_names + model.adsorbate_names, epsilon=1e-100, run_ode=True)
    #model.solver.get_GXRC(17, cvgs, ['C8H6NH-OH_a', 'C8H6NH-H_a', 'C8H6NHOH_a', 'C8H6NH2_a', 'H-H_c'], epsilon=1e-50, run_ode=True)
    #model.solver.get_GXRC(18, cvgs, model.transition_state_names, epsilon=1e-50, run_ode=True)

