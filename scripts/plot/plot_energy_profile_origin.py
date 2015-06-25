import sys
sys.path.append('D:\Dropbox\Code\Python\kinetic\Pynetics')
from pynetics import model
import matplotlib.pyplot as plt

#create micro kinetic model instance
m = model.KineticModel(setup_file='formic_acid.mkm')

energy_tups1 = [(3, 2), (2, 3, 2.5), (1, 4, 3), (3, 4)]
energy_tups2 = [(3, 1.5), (2, 5, 3.5), (3, 4, 2), (3, 4.5)]

energy_tups1 = [m.plotter.get_relative_energy_tuple(energy_tup)
                for energy_tup in energy_tups1]
energy_tups2 = [m.plotter.get_relative_energy_tuple(energy_tup)
                for energy_tup in energy_tups2]

#plot elementary plots
for tups_idx, energy_tups in enumerate((energy_tups1, energy_tups2)):
    for tup_idx, energy_tup in enumerate(energy_tups):
        fname = str(tups_idx)+'_'+str(tup_idx)
        print 'Plotting diagram ' + fname + '...'
        m.plotter.plot_single_energy_diagram(
            rxn_equation=m.rxn_expressions[tup_idx],
            show_mode='save',
            custom_energy_tuple=energy_tup,
            fname=fname
        )
        print 'Ok.'


#plot multi_plots
fig_list = []
points_list = []
for tups_idx, energy_tups in enumerate((energy_tups1, energy_tups2)):
    fname = 'multi_'+str(tups_idx)
    print 'Plotting multi-diagram ' + fname
    fig, x_total, y_total = m.plotter.plot_multi_energy_diagram(
        rxn_equations_list=m.rxn_expressions,
        show_mode='save',
        custom_energy_tuples=energy_tups,
        fname=fname
    )
    fig_list.append(fig)
    points_list.append((x_total, y_total))
    print 'Ok.'

#merge lines
print 'Merge diagrams...'
#lines = []
#for fig in fig_list:
#    lines.extend(fig.axes[0].lines)
new_fig = plt.figure()
ax = new_fig.add_subplot(111)
ax.plot(*points_list[0], linewidth=3, color='#A52A2A')
#m.plotter.add_line_shadow(ax, *points_list[0],
#                          depth=8, color='#595959')
ax.plot(*points_list[1], linewidth=3, color='#00688B')
#m.plotter.add_line_shadow(ax, *points_list[1],
#                          depth=8, color='#595959')
new_fig.show()
print 'Ok.'
