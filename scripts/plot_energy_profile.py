import sys
sys.path.append('D:\Dropbox\Code\Python\kinetic')
from pynetics import model

#create micro kinetic model instance
m = model.KineticModel(setup_file='formic_acid.mkm')

energy_tups1 = []
energy_tups2 = []

energy_tups1 = [m.get_relative_energy_tuple(energy_tup)
                for energy_tup in energy_tups1]
energy_tups2 = [m.get_relative_energy_tuple(energy_tup)
                for energy_tup in energy_tups2]

#plot elementary plots
for tups_idx, energy_tups in enumerate((energy_tups1, energy_tups2)):
    for tup_idx, energy_tup in enumerate(energy_tups):
        m.plotter.plot_single_energy_diagram(
            rxn_equation=m.rxn_expressions[tup_idx],
            show_mode='save',
            custom_energy_tuple=energy_tup,
            fname=str(tups_idx)+'_'+str(tup_idx)
        )

#plot multi_plots
fig_list = []
for tups_idx, energy_tups in enumerate((energy_tups1, energy_tups2)):
    fig = m.plot_multi_energy_diagram(
        rxn_equations_list=m.elementary_rxns_list,
        show_mode='save',
        custom_energy_tuples=energy_tups,
        fname='multi_'+str(tups_idx)
    )
    fig_list.append(fig)

#merge lines
lines = []
for fig in fig_list:
    lines.extend(fig.axes[0].lines)
fig_list[0].axes[0].lines = lines
fig_list[0].show()
