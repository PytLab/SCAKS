import sys
import os
import cPickle as cpkl
import logging

sys.path.append('D:\Dropbox\Code\CentOS_code\Pynetics')
from pynetics import model

if os.path.exists('out.log'):  # del old log file
    print 'delete [ out.log ]...'
    os.remove('out.log')
    print 'Ok.'


if os.path.exists('formulas.tex'):  # del old tex file
    print 'delete [ formulas.tex ]...'
    os.remove('formulas.tex')
    print 'Ok.'

#create micro kinetic model instance
m = model.KineticModel(setup_file='model_setup.mkm')

if len(sys.argv) > 1 and '--correct' in sys.argv:
    correct_energy = True
else:
    correct_energy = False

#use custom initial guess
init_cvgs = (0.0, 0.0, 0.0, 0.0, 0.0)
relative = True if '--relative' in sys.argv else False
solve_ode = True if '--ode' in sys.argv else False
coarse_guess = False if '--nocoarse' in sys.argv else True
m.run_mkm(init_cvgs=init_cvgs, relative=relative,
          correct_energy=correct_energy, solve_ode=solve_ode,
          coarse_guess=True)

# plot energy profiles
#for i, rxn_equation in enumerate(m.rxn_expressions):
#    logging.info('Plotting diagram '+str(i)+'...')
#    m.plotter.plot_single_energy_diagram(rxn_equation, show_mode='save')
#    logging.info('Ok.')
#
#logging.info('Plotting multi_energy_diagram...')
#m.plotter.plot_multi_energy_diagram(m.rxn_expressions, show_mode='save')
#logging.info('Ok.\n')

#get latex file
#logging.info('Generating TEX file...')
#m.solver.get_data_symbols()
#m.solver.get_delta_G_symbols(log_latex=True)
#m.solver.get_rate_syms(log_latex=True)
#m.solver.get_dtheta_dt_syms(log_latex=True)
#logging.info('Ok.')
