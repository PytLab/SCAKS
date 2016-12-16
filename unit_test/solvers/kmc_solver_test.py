import logging
import re
import unittest

import numpy as np

from kynetix.models.kmc_model import KMCModel
from kynetix.solvers import *

from unit_test import *


class KMCSolverTest(unittest.TestCase):

    def setUp(self):
        # Test case setting.
        # {{{
        self.maxDiff = None
        self.setup_dict = dict(
            rxn_expressions = [
                'CO_g + *_t -> CO_t',
                'CO_g + *_b -> CO_b',
                'O2_g + 2*_b -> 2O_b',
                'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b',
                'CO_b + *_t <-> CO_t + *_b -> CO_b + *_t',
            ],

            species_definitions = {
                'CO_g': {'pressure': 0.01},
                'O2_g': {'pressure': 0.2},
                'CO2_g': {'pressure': 0.01},
                'b': {'site_name': 'bridge', 'type': 'site', 'total': 0.5},
                't': {'site_name': 'top', 'type': 'site', 'total': 0.5},
            },

            temperature = 298.,
            parser = "KMCParser",
            solver = "KMCSolver",
            corrector = "ThermodynamicCorrector",
            cell_vectors = [[3.0, 0.0, 0.0],
                            [0.0, 3.0, 0.0],
                            [0.0, 0.0, 3.0]],
            basis_sites = [[0.0, 0.0, 0.0],
                           [0.5, 0.0, 0.0],
                           [0.0, 0.5, 0.0],
                           [0.5, 0.5, 0.0]],
            unitcell_area = 9.0e-20,
            active_ratio = 4./9,
            repetitions = (3, 3, 1),
            periodic = (True, True, False),
            possible_element_types = ["O", "V", "O_s", "C"],
            empty_type = "V",
            possible_site_types = ["P"],
            nstep = 50,
            random_seed = 13996,
            random_generator = 'MT',
            trajectory_dump_interval = 10,
        )
        # }}}

    def test_construction(self):
        " Make sure KMCSolver object can be constructed correctly. "
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        self.assertTrue(isinstance(solver, KMCSolver))

    def test_get_control_parameter(self):
        " Make sure we can get KMCControlParameter object. "
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        solver = model.solver

        control_parameters = solver.get_control_parameters()

        # Check.
        self.assertEqual(10, control_parameters.dumpInterval())
        self.assertEqual(1, control_parameters.analysisInterval())
        self.assertEqual(50, control_parameters.numberOfSteps())
        self.assertEqual(0, control_parameters.rngType())
        self.assertEqual(13996, control_parameters.seed())
        self.assertEqual(False, control_parameters.timeSeed())

    def test_get_rxn_rates_CT(self):
        " Make sure we can get correct forward and reverse rates for a reaction. "
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(relative=True,
                                energy_file=kmc_energy,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)

        ref_r = (1575287.974387463, 3.8789566422291146e-14)
        ret_r = model.solver.get_rxn_rates_CT('CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b', model.relative_energies)
        self.assertTupleEqual(ref_r, ret_r)

        ref_r = (215.85343473385328, 1.7062993852898129e-44)
        ret_r = model.solver.get_rxn_rates_CT('O2_g + 2*_b -> 2O_b', model.relative_energies)
        self.assertTupleEqual(ref_r, ret_r)

        ref_r = (11.535554738754854, 1.3130247359797898e-18)
        ret_r = model.solver.get_rxn_rates_CT('CO_g + *_t -> CO_t', model.relative_energies)
        self.assertTupleEqual(ref_r, ret_r)

    def test_get_single_process(self):
        " Make sure we can parse a process dict correctly. "
        # {{{
        # Construction.
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(relative=True,
                                energy_file=kmc_energy,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)

        process_dict = {"reaction": "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b",
                        "description": "CO and O couple and desorption.",
                        "coordinates_group":[[[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]],
                                             [[0.0, 0.0, 0.0], [0.5, -0.5, 0.0]],
                                             [[0.0, 0.0, 0.0], [-0.5, 0.5, 0.0]],
                                             [[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0]]],
                        "elements_before": ["V", "V"],
                        "elements_after": ["O_s", "O_s"],
                        "basis_sites": [1, 2],
                        "fast": True}

        processes = model.solver._KMCSolver__get_single_process(process_dict)

        # Check processes number.
        self.assertEqual(16, len(processes))

        # Check a the first process object.
        p = processes[0]
        self.assertListEqual(p.basisSites(), [1])
        self.assertListEqual(p.elementsBefore(), ["V", "V"])
        self.assertListEqual(p.elementsAfter(), ["O_s", "O_s"])
        self.assertTrue(p.fast())

        # Check coordinates.
        ref_coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]
        c = p.localConfigurations()[0]
        ret_coords = c.coordinates().tolist()
        self.assertListEqual(ref_coords, ret_coords)

        # Check a the second process object.
        p = processes[-1]
        self.assertListEqual(p.basisSites(), [2])
        self.assertListEqual(p.elementsBefore(), ["O_s", "O_s"])
        self.assertListEqual(p.elementsAfter(), ["V", "V"])
        self.assertTrue(p.fast())

        # Check coordinates.
        ref_coords = [[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0]]
        c = p.localConfigurations()[0]
        ret_coords = c.coordinates().tolist()
        self.assertListEqual(ref_coords, ret_coords)
        # }}}

    def test_run(self):
        " Test the we can run the kmc model correctly. "
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(energy_file=kmc_energy,
                                relative=True,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)
        solver = model.solver

        solver.run()

    def test_process_reactions_mapping(self):
        " Test process mapping query function. "
        # {{{
        model = KMCModel(setup_dict=self.setup_dict, logger_level=logging.WARNING)
        model.parser.parse_data(energy_file=kmc_energy,
                                relative=True,
                                processes_file=kmc_processes,
                                configuration_file=kmc_config,
                                sitesmap_file=kmc_sites)
        processes = model.solver.processes

        for p in processes:
            self.assertEqual(p.__class__.__name__, "KMCProcess")

        ret_mapping = model.solver.process_mapping

        ref_mapping = ['CO_g + *_t -> CO_t(->)',
                       'CO_g + *_t -> CO_t(<-)',
                       'CO_g + *_b -> CO_b(->)',
                       'CO_g + *_b -> CO_b(<-)',
                       'CO_g + *_b -> CO_b(->)',
                       'CO_g + *_b -> CO_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'O2_g + 2*_b -> 2O_b(->)',
                       'O2_g + 2*_b -> 2O_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(->)',
                       'CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b(<-)',]

        self.assertListEqual(ref_mapping, ret_mapping)
        # }}}

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(KMCSolverTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

