import sys
import time

def fmtcovert(sec):
    hours = int(sec/(3600.0))
    minutes = int((sec - hours*3600)/60.0)
    seconds = sec - hours*3600 - minutes*60

    return hours, minutes, seconds  # int, int, float

sys.path.append('D:\Dropbox\Code\CentOS_code\Pynetics')
from pynetics import model
if '--traj' not in sys.argv:
    from pynetics.solvers.kmc_plugins import *

#create micro kinetic model instance
m = model.KineticModel(setup_file='kmc_setup.mkm')

if '--traj' not in sys.argv:
    start = time.time()
    m.solver.run()
    end = time.time()

    t = end - start
    h, m, s = fmtcovert(t)
    print "Time used: %d h %d min %f sec" % (h, m, s)
