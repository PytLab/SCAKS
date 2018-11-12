import inspect
import os
from ..compatutil import subprocess

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
mkm_abs_energy = mkm_path + "/abs_energy.py"

def cleanup():
    ''' Remove auto-generated files.
    '''
    subprocess.getstatusoutput("for i in `find ./ -name 'auto_*'`; do rm -rf $i; done")
    subprocess.getstatusoutput("for i in `find ./ -name 'out.log'`; do rm -rf $i; done")
    subprocess.getstatusoutput("for i in `find ./ -name 'log'`; do rm -rf $i; done")
    subprocess.getstatusoutput("for i in `find ./ -name '*.pkl'`; do rm -rf $i; done")
    subprocess.getstatusoutput("rm ./*.csv")

