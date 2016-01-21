import sys
import os
import cPickle as cpkl

sys.path.append('D:\Dropbox\Code\CentOS_code\Pynetics')
from pynetics import model


#create micro kinetic model instance
m = model.KineticModel(setup_file='model_setup.mkm')

if len(sys.argv) > 1:
    if sys.argv[1] == '--init':  # table
        m.table_maker.create_initial_table()
    else:
        if '--update' in sys.argv:  # update
            m.table_maker.create_new_table('update')
        if '--add' in sys.argv:  # add
            m.table_maker.create_new_table('add')

        m.parser.parse_data()  # parse data from table
        m.solver.get_data_dict()  # solver get data from parser
        m.solver.get_rate_constants()
        m.solver.get_rate_expressions(m.solver.rxns_list)

        for i, rxn_equation in enumerate(m.rxn_expressions):
            print 'Plotting diagram '+str(i)+'...'
            m.plotter.plot_single_energy_diagram(rxn_equation, show_mode='save')
            print 'Ok.'

        print "Plotting multi_energy_diagram..."
        m.plotter.plot_multi_energy_diagram(m.rxn_expressions, show_mode='save')
        print 'Ok.\n'

        if 'c' in sys.argv[1]:
            print 'Correct free energies...'
            m.solver.correct_energies()
            print 'Ok.'

        #set initial guess(initial coverage)
        #init_cvg = (0.0, 0.0, 0.0, 0.0, 0.0, 0.00, 0.0, 0.0, 0.99, 0.0, 0.0, 0.0, 0.0, 0.0)
        #if there is converged coverage, use it as initial guess
        if os.path.exists("./data.pkl"):
            with open('data.pkl', 'rb') as f:
                data = cpkl.load(f)
            if 'steady_state_coverage' in data:
                init_cvg = data['steady_state_coverage']
            elif "initial_guess" in data:
                init_cvg = data['initial_guess']
            else:
                init_cvg = m.solver.boltzmann_coverages()
        else:  # use Boltzmann coverage
            init_cvg = m.solver.boltzmann_coverages()
        cvg = m.solver.get_steady_state_cvgs(init_cvg)

        #get latex file
        print "Generating TEX file..."
        m.solver.get_data_symbols()
        m.solver.get_delta_G_symbols(log_latex=True)
        m.solver.get_rate_syms(log_latex=True)
        m.solver.get_dtheta_dt_syms(log_latex=True)
        print "OK."
