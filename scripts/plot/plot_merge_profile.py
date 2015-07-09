'''
    Script to plot merged energy profile
'''
import sys

from simple_plot import *
from merge_data import *  # get rxn_equations & energy_tuples
import matplotlib.pyplot as plt

points_list = []
print "Plotting single multi-energy profile..."
for idx, energy_tuples in enumerate(multi_energy_tuples):
    fname = 'multi_energy_diagram_' + str(idx).zfill(2)
    print "Plotting diagram " + fname + "..."
    fig, x_total, y_total = \
        plot_multi_energy_diagram(rxn_equations, energy_tuples,
                                  n=10000, show_mode='save', fname=fname)
    print "Ok."
    points_list.append((x_total, y_total))

#merge lines
print 'Merge diagrams...'
new_fig = plt.figure(figsize=(16, 9))
# transparent figure
if sys.argv[2] == '--trans':
    new_fig.patch.set_alpha(0)

ax = new_fig.add_subplot(111)
# transparent axe
if sys.argv[2] == '--trans':
    ax.patch.set_alpha(0)

#remove xticks
ax.set_xticks([])
ax.set_xmargin(0.03)
ax.set_ylim(-2.25, 0.25)
ax.set_yticks(np.linspace(-2.25, 0.25, 11))
ax.set_yticklabels(['', '-2.0', '', '-1.5', '', '-1.0', '', '-0.5', '', '0.0', ''])

for color, points in zip(colors, points_list):
    add_line_shadow(ax, *points, depth=7, color='#595959',
                    line_width=3, offset_coeff=1.0)
    ax.plot(*points, linewidth=4.5, color=color)
if sys.argv[1] == '--show':
    new_fig.show()
elif sys.argv[1] == '--save':
    new_fig.savefig('merged_energy_profile.png', dpi=1000)
print 'Ok.'
