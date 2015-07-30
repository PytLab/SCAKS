'''
    Script to plot energy profile
'''
import sys

from simple_plot import *
from data import *  # get rxn_equations & energy_tuples

if 'rxn_equations' in dir() and 'energy_tuples' in dir():  # multi rxn
    #check data shape
    if len(rxn_equations) != len(energy_tuples):
        raise ValueError("lengths of rxn_equations and energy_tuples " +
                         "are different.")
    for rxn_equation, energy_tuple in zip(rxn_equations, energy_tuples):
        equation_list = equation2list(rxn_equation)
        if len(equation_list) != len(energy_tuple):
            raise ValueError("unmatched shape: %d, %d" %
                             (rxn_equation, str(energy_tuple)))
    #plot single diagrams
    for idx, args in enumerate(zip(energy_tuples, rxn_equations)):
        fname = str(idx).zfill(2)
        print "Plotting diagram " + fname + "..."
        plot_single_energy_diagram(*args, show_mode='save', fname=fname)
        print "Ok."
    #plot multi-diagram
    print "Plotting multi-diagram..."
    fig, x_total, y_total = \
        plot_multi_energy_diagram(rxn_equations, energy_tuples, show_mode='save')
    print "Ok."
elif 'rxn_equation' in dir() and 'energy_tuple' in dir():  # single rxn
    print "Plotting multi-diagram..."
    fig, x_total, y_total = \
        plot_single_energy_diagram(energy_tuple, rxn_equation, show_mode='save')
    print "Ok."
else:  # no equation and energy tuple
    raise ValueError('No rxn equation and energy tuple is defined.\n' +
                     'Please check you data file...')

#customize your diagram
if 'custom' in dir() and custom:
    print "Custom plotting..."
    new_fig = plt.figure(figsize=(16, 9))
    # transparent figure
    if len(sys.argv) > 2 and sys.argv[2] == '--trans':
        new_fig.patch.set_alpha(0)

    ax = new_fig.add_subplot(111)
    # transparent axe
    if len(sys.argv) > 2 and sys.argv[2] == '--trans':
        ax.patch.set_alpha(0)
    #remove xticks
    ax.set_xticks([])
    ax.set_xmargin(0.03)

    #set attributes of y-axis
    if 'ylim' in dir():
        ymin, ymax = ylim
        ax.set_ylim(ymin, ymax)
        if 'n_yticks' in dir():  # must set ylim befor setting n_yticks
            ax.set_yticks(np.linspace(ymin, ymax, n_yticks))
    if 'yticklabels' in dir():
        ax.set_yticklabels(yticklabels)

    #add line shadow
    shadow_depth = shadow_depth if 'shadow_depth' in dir() else 7
    shadow_color = shadow_color if 'shadow_color' in dir() else '#595959'
    offset_coeff = offset_coeff if 'offset_coeff' in dir() else 9.0
    add_line_shadow(ax, x_total, y_total, depth=shadow_depth,
                    color=shadow_color, line_width=5.4,
                    offset_coeff=offset_coeff)

    if 'color' not in dir():
        print "No custom color. \nUse default color: black."
        color = '#000000'
    ax.plot(x_total, y_total, linewidth=5.4, color=color)
    if sys.argv[1] == '--show':
        new_fig.show()
    elif sys.argv[1] == '--save':
        new_fig.savefig('./energy_profile/energy_profile.png', dpi=500)
    else:
        raise ValueError('Unrecognized show mode parameter : %s.', sys.argv[1])
    print 'Ok.'
