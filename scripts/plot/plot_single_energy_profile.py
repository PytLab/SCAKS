'''Script to plot energy profile'''
from simple_plot import *
from single_data import *  # get rxn_equations & energy_tuples

#check data shape
if len(rxn_equations) != len(energy_tuples):
    raise ValueError("lengths of rxn_equations and energy_tuples " +
                     "are different.")

for idx, args in enumerate(zip(energy_tuples, rxn_equations)):
    fname = str(idx).zfill(2)
    print "Plotting diagram " + fname + "..."
    plot_single_energy_diagram(*args, show_mode='save', fname=fname)
    print "Ok."
