import numpy as np

from scaks.models.micro_kinetic_model import MicroKineticModel
from scaks.solvers.steady_state_solver import SteadyStateSolver
from scaks.correctors.thermodynamic_corrector import ThermodynamicCorrector

model_dict = {
    'rxn_expressions': ['CO_g + *_s -> CO_s',
                        'O2_g + 2*_s <-> O-O_2s -> 2O_s',
                        'CO_s + O_s <-> CO-O_2s -> CO2_g + 2*_s'],

    'species_definitions': {
        'CO_g': {'pressure': 1.32e-7},
        'O2_g': {'pressure': 5.26e-7},
        'CO2_g': {'pressure': 1.32e-7},
        '*_s': {'site_name': 'Pt100', 'type': 'site', 'total': 1.0},
    },

    'temperature': 500,
    'parser': "RelativeEnergyParser",

    'rate_algo': 'CT',
    'unitcell_area': 9.0e-20,
    'active_ratio': 4.0/9.0,

    'tolerance': 1e-50
}

relative_energies = {
    'Ga': [0.0, 0.07, 0.39],
    'dG': [-2.09, -2.39, -0.46]
}

# Build micor-kinetic model.
model = MicroKineticModel(setup_dict=model_dict)

# Read data.
model.parser.parse_data(energy_data=relative_energies)

# Create solver
solver = SteadyStateSolver(model)
model.set_solver(solver)
model.solver.get_data()

# Create corrector for energy correction
corrector = ThermodynamicCorrector(model)
model.set_corrector(corrector)

# Create and register a custom hybrid method
@model.hybrid_method_register
def random_guess(model, N):
    if model.log_allowed:
        model.logger.info('Re-initialize coverages randomly')

    n_adsorbates = len(model.adsorbate_names)
    new_cvgs = np.random.randn(n_adsorbates).tolist()

    return new_cvgs

if __name__ == '__main__':

    # Run.
    model.run(init_cvgs=[0.0, 0.0],
              solve_ode=False,
              coarse_guess=False,
              XRC=False,
              product_name='CO2_g')

