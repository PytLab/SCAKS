import inspect
import os

abs_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

# Get absolute path for kmc input data.
kmc_path = abs_path + "/kmc_inputs"
kmc_energy = kmc_path + "/rel_energy.py"
kmc_processes = kmc_path + "/kmc_processes.py"
kmc_sites = kmc_path + "/kmc_sites.py"
kmc_config = kmc_path + "/kmc_configuration.py"

# Get absolute path for mkm input data.
mkm_path = abs_path + "/mkm_inputs"
mkm_energy = mkm_path + "/rel_energy.py"

