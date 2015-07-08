'''
    Script to plot energy profile
'''
from simple_plot import *
from data import *  # get rxn_equations & energy_tuples

if 'rxn_equations' in dir() and 'energy_tuples' in dir():  # single rxn or multi rxn
    #check data shape
    if len(rxn_equations) != len(energy_tuples):
        raise ValueError("lengths of rxn_equations and energy_tuples " +
                         "are different.")

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
else:
    print "Plotting multi-diagram..."
    fig, x_total, y_total = \
        plot_single_energy_diagram(energy_tuple, rxn_equation, show_mode='save')
    print "Ok."

#if u want to customize your diagram
#remove the comment symbols and modify codes below
'''
print "Custom plotting..."
new_fig = plt.figure(figsize=(16, 9))
ax = new_fig.add_subplot(111)
#remove xticks
ax.set_xticks([])
ax.set_ylim(-0.75, 0.5)
ax.set_yticks(np.linspace(-0.75, 0.5, 11))
ax.set_yticklabels(['', '', '-0.5', '', '', '', '0.0', '', '', '', '0.5'])
ax.set_xmargin(0.03)
#add line shadow
add_line_shadow(ax, x_total, y_total, depth=7, color='#595959', line_width=6, offset_coeff=12.0)
ax.plot(x_total, y_total, linewidth=6, color='#A52A2A')
new_fig.savefig('energy_profile.png', dpi=1200)
#new_fig.show()
print 'Ok.'
'''
