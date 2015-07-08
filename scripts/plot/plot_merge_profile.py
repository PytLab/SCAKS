'''
    Script to plot merged energy profile
'''
from simple_plot import *
from merge_data import *  # get rxn_equations & energy_tuples
import matplotlib.pyplot as plt

points_list = []
for energy_tuples in multi_energy_tuples:
    fig, x_total, y_total = \
        plot_multi_energy_diagram(rxn_equations, energy_tuples,
                                  n=10000, show_mode='save')
    points_list.append((x_total, y_total))

#merge lines
print 'Merge diagrams...'
new_fig = plt.figure(figsize=(16, 9))
ax = new_fig.add_subplot(111)
#remove xticks
ax.set_xticks([])
ax.set_xmargin(0.03)
ax.set_ylim(-2.25, 0.25)
ax.set_yticks(np.linspace(-2.25, 0.25, 11))
ax.set_yticklabels(['', '-2.0', '', '-1.5', '', '-1.0', '', '-0.5', '', '0.0', ''])
add_line_shadow(ax, *points_list[0], depth=7, color='#595959', line_width=3, offset_coeff=12.0)
ax.plot(*points_list[0], linewidth=4.5, color='#A52A2A')
add_line_shadow(ax, *points_list[1], depth=7, color='#595959', line_width=3, offset_coeff=12.0)
ax.plot(*points_list[1], linewidth=4.5, color='#000000')
#new_fig.show()
new_fig.savefig('merged_energy_profile.png', dpi=1000)
print 'Ok.'
