import sys
import os
import cPickle as cpkl

sys.path.append('D:\Dropbox\Code\Python\kinetic\Pynetics')
from pynetics import model

if os.path.exists('event.log'):  # del old log file
    os.remove('event.log')
    print "delete 'event.log'."
#create micro kinetic model instance
m = model.KineticModel(setup_file='furfural.mkm')

m.parser.parse_data()  # parse data from rel_energy.py
m.solver.get_data_dict()  # solver get data from parser
m.solver.get_rate_constants()
m.solver.get_rate_expressions(m.solver.rxns_list)
'''
for i, rxn_equation in enumerate(m.rxn_expressions):
    print 'Plotting diagram '+str(i)+'...'
    m.plotter.plot_single_energy_diagram(rxn_equation, show_mode='save')
    print 'Ok.'

print "Plotting multi_energy_diagram..."
m.plotter.plot_multi_energy_diagram(m.rxn_expressions, show_mode='save')
print 'Ok.\n'
'''

if len(sys.argv) > 1 and 'c' in sys.argv[1]:
    print 'Correct free energies...'
    m.solver.correct_energies()
    print 'Ok.'

#set initial guess(initial coverage)
#init_cvg = (0.0, 0.0, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.99, 0.0, 0.0, 0.0, 0.0, 0.0)
#if there is converged coverage, use it as initial guess
if os.path.exists("./data.pkl"):
    with open('data.pkl', 'rb') as f:
        data = cpkl.load(f)
    #init_guess = 'initial_guess'
    init_guess = 'steady_state_coverage'
    if init_guess in data:
        init_cvg = data[init_guess]
    else:
        init_cvg = m.solver.boltzmann_coverages()
else:  # use Boltzmann coverage
    init_cvg = m.solver.boltzmann_coverages()

init_cvg = (0.0, 0.0, 0.99, 0.0, 0.0, 0.00, 0.0, 0.0, 0.00, 0.0, 0.0, 0.0, 0.0, 0.0)

cvg = m.solver.get_steady_state_cvgs(init_cvg)

#get latex file
print "Generating TEX file..."
m.solver.get_data_symbols()
m.solver.get_delta_G_symbols(log_latex=True)
m.solver.get_rate_syms(log_latex=True)
m.solver.get_dtheta_dt_syms(log_latex=True)
print "OK."
