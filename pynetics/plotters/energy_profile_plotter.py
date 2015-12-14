import threading
import logging

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from math import sqrt
from matplotlib.patches import Ellipse

from plotter_base import PlotterBase


class ShadowThread(threading.Thread):
    "Sub-class of Thread class to create threads to plot shadows."
    def __init__(self, func, args):
        threading.Thread.__init__(self)
        self.func = func
        self.args = args

    def run(self):
        apply(self.func, self.args)


class NoteThread(threading.Thread):
    "Sub-class of Thread class to add notes to line."
    def __init__(self, func, args):
        threading.Thread.__init__(self)
        self.func = func
        self.args = args

    def run(self):
        apply(self.func, self.args)


class EnergyProfilePlotter(PlotterBase):
    def __init__(self, owner):
        PlotterBase.__init__(self, owner)
        # set logger
        self.logger = logging.getLogger('model.plotters.EnProfilePlotter')

    @staticmethod
    def quadratic_interp_poly(x1, y1, x2, y2):
        """
        Return polynomial function object f(x)
        obtained from quadratic interpolation.

        Parameters
        ----------
        x1, y1 : float
            x and y value of the left point.
        x2, y2 : float
            x and y value of the left point(maximum value point).

        Examples
        --------
        >>> f = m.plotter.quadratic_interp_poly(0.0, 0.0, 2.0, 2.0)
        >>> -0.5, 2, 0, <function kinetic.plotters.general_plotter.poly_func>

        """
        A = np.matrix([[x1**2, x1, 1],
                       [x2**2, x2, 1],
                       [ 2*x2,  1, 0]])

        b = np.matrix([[y1], [y2], [0]])
        x = A.I * b
        x.shape = (1, -1)
        a, b, c = x.tolist()[0]

        def poly_func(x):
            return a*x**2 + b*x + c

        return a, b, c, poly_func

    def get_potential_energy_points(self, energy_tuple, n=100,
                                    subsection_length=1.0):
        """
        Expect a energuy_tuple containing 3 components, number of points
        return two 1-D arrays(x, y).

        Parameters
        ----------
        energy_tuple : tuple
            A tuple containing (E_IS, E_TS, E_FS).
        n : float, optional
            Number of points interpolated.
        subsection_length : int, optional
            Length of each subsection(x_i, x_f)
        """
        #use quadratic interpolation to get barrier points
        if len(energy_tuple) == 3:
            #transition state y
            a, b, c, f = self.quadratic_interp_poly(0.0, energy_tuple[0],
                                                    subsection_length/2,
                                                    energy_tuple[1])
            y3 = energy_tuple[-1]
            try:
                x3 = max((-b + sqrt(b**2 - 4*a*(c - y3)))/(2*a),
                         (-b - sqrt(b**2 - 4*a*(c - y3)))/(2*a))
            except ValueError:
                raise ValueError('abnormal energy : ' + str(energy_tuple))
            init_x_b = np.linspace(0, x3, n)
            f_ufunc = np.frompyfunc(f, 1, 1)  # convert to universal function
            y_b = f_ufunc(init_x_b)
            x_b = np.linspace(subsection_length, subsection_length+x3, n)
            #initial state y
            x_i = np.linspace(0, subsection_length, n)
            y_i = np.linspace(energy_tuple[0], energy_tuple[0], n)
            #final state y
            x_f = np.linspace(subsection_length+x3, 2*subsection_length+x3, n)
            y_f = np.linspace(energy_tuple[-1], energy_tuple[-1], n)
            #combine all points
            y = np.array(y_i.tolist() + y_b.tolist() + y_f.tolist())
            #x = np.linspace(0, x3 + 2*subsection_length, 3*n)
            x = np.array(x_i.tolist() + x_b.tolist() + x_f.tolist())
            ts_scale = x3

        if len(energy_tuple) == 2:
            #transition state
            energy_list = list(energy_tuple)
            if energy_list[0] < energy_list[-1]:
                energy_list.insert(1, energy_tuple[-1]+1e-100)
                init_x_b = np.array([0.0, subsection_length - 1e-5,
                                     subsection_length])
            else:
                energy_list.insert(1, energy_tuple[0]+1e-100)
                init_x_b = np.array([0.0, 1e-5, subsection_length])
            energy_tuple = tuple(energy_list)

            init_y_b = np.array(energy_tuple)
            f_b = interpolate.interp1d(init_x_b, init_y_b, kind='quadratic')
            x_b = np.linspace(0, subsection_length, n)
            y_b = f_b(x_b)
            #initial state y
            y_i = np.linspace(energy_tuple[0], energy_tuple[0], n)
            #final state y
            y_f = np.linspace(energy_tuple[-1], energy_tuple[-1], n)
            #combine all points
            y = np.array(y_i.tolist() + y_b.tolist() + y_f.tolist())
            x = np.linspace(0, 3*subsection_length, 3*n)
            ts_scale = subsection_length

        #plt.plot(x, y)
        #plt.show()
        return x, y, ts_scale

    def get_energy_tuple(self, rxn_equation):
        """
        Expect a reaction equation,
        return its energy tuple which contains(E_IS, E_TS, E_FS) or (E_IS, E_FS),
        set E_IS as 0.0 and other energies are relative energy to E_IS.
        """
        idx = self._owner.rxn_expressions.index(rxn_equation)
        elementary_rxn_list = self._owner.elementary_rxns_list[idx]
        if not hasattr(self._owner, 'state_energy_dict'):
            getattr(self._owner.solver, 'get_state_energy_dict')()
        #go through elementary rxn list to get energy tuple
        energy_list = []
        for state in elementary_rxn_list:
            state_str = ' + '.join(state)
            energy_list.append(float(self._owner.state_energy_dict[state_str]))
        #set initial energy to 0 (reference)
        reference_energy = energy_list[0]
        for i in xrange(len(energy_list)):
            energy_list[i] = energy_list[i] - reference_energy
        #compare E_IS and E_TS
        if len(energy_list) == 3:
            if energy_list[1] <= energy_list[0] or \
               energy_list[1] <= energy_list[-1]:
                init_energy_tup = tuple(energy_list)
                energy_list = [energy_list[0], energy_list[-1]]
                #log the event
                self.logger.warning('Unreasonable energy tuple %s',
                                    str(init_energy_tup))
                self.logger.warning('               Convert to %s',
                                    str(tuple(energy_list)))

        return tuple(energy_list)

    @staticmethod
    def get_relative_energy_tuple(energy_tuple):
        energy_list = list(energy_tuple)
        reference_energy = energy_list[0]
        for i in xrange(len(energy_list)):
            energy_list[i] = energy_list[i] - reference_energy
        return tuple(energy_list)

    @staticmethod
    def add_line_shadow(ax, x, y, depth, color, line_width=3):
        "add shadow to line in axes 'ax' by changing attribute of the object"
        def add_single_shadow(ax, x, y, order, depth, color, line_width):
            offset = transforms.ScaledTranslation(order, -order,
                                                  transforms.IdentityTransform())
            shadow_trans = ax.transData + offset
            ax.plot(x, y, linewidth=line_width, color=color,
                    transform=shadow_trans,
                    alpha=(depth-order)/2.0/depth)
            return

        threads = []
        for i in xrange(depth):  # gather thread objects
            t = ShadowThread(add_single_shadow, (ax, x, y, i, depth,
                                                 color, line_width))
            threads.append(t)

        for i in xrange(depth):
            threads[i].start()

        for i in xrange(depth):
            threads[i].join()

        return ax.lines

    def plot_single_energy_diagram(self, rxn_equation, n=100,
                                   subsection_length=1.0,
                                   line_color='#000000',
                                   has_shadow=True, fmt='jpeg',
                                   show_mode='show', **kwargs):
        """
        Draw a potential energy diagram of a elementary reaction equation
        and save it.
        Return the Figure object.

        Parameters
        ----------
        rxn_equation : str
            reaction equation string according to rules in setup file.
            e.g. 'HCOOH_g + 2*_s <-> HCOO-H_s + *_s -> HCOO_s + H_s'.

        subsection_length : int, optional
            Length of each subsection(x_i, x_f), defualt to be 1.0.

        line_color : str, optional
            Color code of the line. Default to be '#000000'/'black'.

        has_shadow : Bool
            Whether to add shadow effect to the main line.
            Default to be False.

        fmt : str, optional
            The output format : ['pdf'|'jpeg'|'png'|'raw'|'pgf'|'ps'|'svg']
            Default to be 'jpeg'

        show_mode : str, 'show' or 'save'
            'show' : show the figure in a interactive way.
            'save' : save the figure automatically.

        Other parameters
        ----------------
        kwargs :
            'fname' : str, optional
            'custom_energy_tuple' : tuple, optional

        Examples
        --------
        >>> m.plotter.plot_single_energy_diagram('COO_s -> CO2_g + *_s',
                                                 has_shadow=True,
                                                 'fname'='pytlab')
        >>> <matplotlib.figure.Figure at 0x5659f30>
        """
        idx = self._owner.rxn_expressions.index(rxn_equation)
        elementary_rxn_list = self._owner.elementary_rxns_list[idx]
        #get energy tuple
        if 'custom_energy_tuple' in kwargs:
            energy_tuple = kwargs['custom_energy_tuple']
        else:
            energy_tuple = self.get_energy_tuple(rxn_equation)
        #energy info
        rxn_energy = round(energy_tuple[-1] - energy_tuple[0], 2)
        if len(energy_tuple):
            act_energy = round(energy_tuple[1] - energy_tuple[0], 2)
        #get x, y array
        x, y, ts_scale = \
            self.get_potential_energy_points(energy_tuple, n=n,
                                             subsection_length=subsection_length)
        #get maximum and minimum values of x and y axis
        y_min, y_max = np.min(y), np.max(y)
        x_max = 2.7*subsection_length + ts_scale
        #scale of y
        y_scale = abs(y_max - y_min)

        ###################      start Artist Mode!      #######################
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111)
        ax.set_ylim(y_min - y_scale/5, y_max + y_scale/5)
        ax.set_xlim(-subsection_length/3, x_max)
        xlabel = ax.set_xlabel('Reaction Coordinate')
        ylabel = ax.set_ylabel('Free Energy(eV)')
        latex_title = self.get_equation_latex_expression(elementary_rxn_list)
        title = ax.set_title(latex_title,
                             fontdict={
                                 'fontsize': 13,
                                 'weight': 1000,
                                 'verticalalignment': 'bottom'
                                 })

        #add shadow
        if has_shadow:
            self.add_line_shadow(ax, x, y, depth=8, color='#595959')
        #main line
        ax.plot(x, y, color=line_color, linewidth=3)

        #energy latex string
        if act_energy:
            act_energy_latex = r'$\bf{G_{a} = ' + str(act_energy) + r' eV}$'
        rxn_energy_latex = r'$\bf{\Delta G = ' + str(rxn_energy) + r' eV}$'

        #######################   annotates   #######################

        #add species annotation
        #get annotate coordiantes
        note_offset = y_scale/40
        param_list = []
        #initial state
        note_x_i = subsection_length/10
        note_y_i = energy_tuple[0] + note_offset
        note_str_i = r'$\bf{'+self.state_latex_convert(elementary_rxn_list[0])+r'}$'
        param_list.append((note_x_i, note_y_i, note_str_i))
        #final state
        note_x_f = subsection_length/10 + subsection_length + ts_scale
        note_y_f = energy_tuple[-1] + note_offset
        note_str_f = r'$\bf{'+self.state_latex_convert(elementary_rxn_list[-1])+r'}$'
        param_list.append((note_x_f, note_y_f, note_str_f))
        if len(energy_tuple) == 3:
            #transition state
            note_y_b = np.max(y) + note_offset  # equal to energy_tuple[1]
            idx = np.argmax(y)
            note_x_b = x[idx] - ts_scale/3
            #####################
            barrier_x = x[idx]  # x value of TS point
            #####################
            note_str_b = r'$\bf'+self.state_latex_convert(elementary_rxn_list[1])+r'}$'
            param_list.append((note_x_b, note_y_b, note_str_b))
        #add annotates
        for idx, param_tuple in enumerate(param_list):
            if idx == 2:
                ax.text(*param_tuple, fontdict={'fontsize': 13,
                                                'color': '#CD5555'})
            else:
                ax.text(*param_tuple, fontdict={'fontsize': 13,
                                                'color': '#1874CD'})

        ######################################################
        #                                                    #
        #  inflection points coordinates:                    #
        #    (subsection_length, energy_tuple[0])            #
        #    (barrier_x, np.max(y))                          #
        #    (subsection_length+ts_scale, energy_tuple[-1])  #
        #                                                    #
        ######################################################

        #energy change annotation

        #initial state horizontal line
        hori_x_i = np.linspace(subsection_length,
                               2*subsection_length+ts_scale, 50)
        hori_y_i = np.linspace(energy_tuple[0], energy_tuple[0], 50)
        ax.plot(hori_x_i, hori_y_i, color=line_color, linewidth=1,
                linestyle='dashed')

        el = Ellipse((2, -1), 0.5, 0.5)

        if len(energy_tuple) == 3:  # for reaction with barrier
            #arrow between barrier
            ax.annotate('', xy=(barrier_x, energy_tuple[0]), xycoords='data',
                        xytext=(barrier_x, np.max(y)), textcoords='data',
                        arrowprops=dict(arrowstyle="<->"))
            #text annotations of barrier
            ax.annotate(act_energy_latex,
                        xy=(barrier_x, np.max(y)/2), xycoords='data',
                        xytext=(-150, 30), textcoords='offset points',
                        size=13, color='#FF4500',
                        arrowprops=dict(arrowstyle="simple",
                                        fc="0.6", ec="none",
                                        patchB=el,
                                        connectionstyle="arc3,rad=0.2"),
                        )

        #arrow between reaction energy
        ax.annotate('',
                    xy=(subsection_length*3/2+ts_scale, energy_tuple[-1]),
                    xycoords='data',
                    xytext=(subsection_length*3/2+ts_scale, energy_tuple[0]),
                    textcoords='data',
                    arrowprops=dict(arrowstyle="<->")
                    )
        #text annotations of reaction energy
        ax.annotate(rxn_energy_latex, xy=(subsection_length*3/2+ts_scale,
                                          (energy_tuple[-1]+energy_tuple[0])/2),
                    xycoords='data',
                    xytext=(50, 30), textcoords='offset points',
                    size=13, color='#8E388E',
                    arrowprops=dict(arrowstyle="simple",
                                    fc="0.6", ec="none",
                                    patchB=el,
                                    connectionstyle="arc3,rad=-0.2"),
                    )
              ##############      annotates end      ################
         #################      Artist Mode End !      #######################
        #remove xticks
        ax.set_xticks([])
        #save the figure object
        if 'fname' in kwargs:
            fname = kwargs['fname'] + '.' + fmt
        else:
            fname = str(self._owner.rxn_expressions.index(rxn_equation))+'.'+fmt

        if show_mode == 'save':
            fig.savefig(fname)
        if show_mode == 'show':
            plt.show()

        return fig

    def pure_species_latex_convert(self, pure_species_str):
        "e.g. convert 'CO2' to 'CO_{2}'"
        pure_species_list = list(pure_species_str)
        for idx, element in enumerate(pure_species_list):
            try:
                num = float(element)
                pure_species_list[idx] = '_{' + element + '}'
            except ValueError:
                pass
        return ''.join(pure_species_list)

    def species_latex_convert(self, species_str):
        """
        Expect a species_str e.g. 'HCOO_s',
        return a latex string 'HCOO^*(s)'
        """
        stoichiometry, species_name = self.split_species(species_str)
        #for free site
        if '*' in species_name:
            tex_sp = '*(s)'
        else:
            sp, site = species_name.split('_')
            if site == 'g':
                tex_site = '(g)'
            else:
                tex_site = '^*(' + site + ')'
            tex_sp = self.pure_species_latex_convert(sp) + tex_site
        if stoichiometry != 1:
            tex_species = str(stoichiometry) + tex_sp
        else:
            tex_species = tex_sp

        return tex_species

    def state_latex_convert(self, state_list):
        """
        Expect a state_list e.g. ['HCOOH_g', '2*_s'],
        return latex string 'HCOOH(g) + 2*(s)'
        """
        tex_state_list = []
        for species_str in state_list:
            tex_state_list.append(self.species_latex_convert(species_str))
        tex_state = ' +  '.join(tex_state_list)
        return tex_state

    def get_equation_latex_expression(self, elementary_rxn_list):
        tex_states_list = []
        for state_list in elementary_rxn_list:
            tex_state = self.state_latex_convert(state_list)
            tex_states_list.append(tex_state)
        #combine tex_states_list
        state_num = len(tex_states_list)
        if state_num == 3:
            tex_equation = tex_states_list[0] + r' \leftrightarrow  ' + \
                           tex_states_list[1] + r' \rightarrow  ' + \
                           tex_states_list[2]
        if state_num == 2:
            tex_equation = tex_states_list[0] + r' \rightarrow  ' + \
                           tex_states_list[1]

        return r'$\bf{' + tex_equation + r'}$'

    def rxn_string2list(self, rxn_equation):
        idx = self._owner.rxn_expressions.index(rxn_equation)
        return self._owner.elementary_rxns_list[idx]

    def plot_multi_energy_diagram(self, rxn_equations_list, n=100,
                                  subsection_length=1.0, line_color='#000000',
                                  fmt='jpeg', has_shadow=True,
                                  show_note=True, show_aux_line=True,
                                  show_mode='show', show_arrow=True,
                                  **kwargs):
        """
        Draw a potential energy diagram of a series of
        elementary reaction equation and save or show it.
        Return the Figure object.

        Parameters
        ----------
        rxn_equations_list : list
            a list containing a series of reaction equation string
            according to rules in setup file.
            e.g. ['HCOOH_g + 2*_s <-> HCOO-H_s + *_s -> HCOO_s + H_s',
                  'HCOO_s + *_s <-> H-COO_s + *_s -> COO_s + H_s',
                  'COO_s -> CO2_g + *_s',
                  '2H_s <-> H-H_s + *_s -> H2_g + 2*_s'].

        subsection_length : int, optional
            Length of each subsection(x_i, x_f), defualt to be 1.0.

        line_color : str, optional
            Color code of the line. Default to be '#000000'/'black'.

        has_shadow : Bool
            Whether to add shadow effect to the main line.
            Default to be True.

        show_note : bool
            Whether to show species note on the main line.
            Default to be True.

        show_aux_line: bool
            Whether to show auxiliary lines on the main line.
            Default to be Ture.

        show_arrow : bool
            Whether to show annotation arrow.
            Default to be True.

        fmt : str, optional
            The output format : ['pdf'|'jpeg'|'png'|'raw'|'pgf'|'ps'|'svg']
            Default to be 'jpeg'

        show_mode : str, 'show' or 'save'
            'show' : show the figure in a interactive way.
            'save' : save the figure automatically.

        Other parameters
        ----------------
        kwargs :
            'fname' : str, optional
            'custom_energy_tuples': list of tuples

        Examples
        --------

        """
        #convert rxn_equation_list(str) to elementary_rxns_list alike list
        rxns_list = [self.rxn_string2list(rxn_equation)
                     for rxn_equation in rxn_equations_list]

        x_offset, y_offset = 0.0, 0.0
        total_x, total_y = [], []
        piece_points, ts_scales = [], []  # collectors

        for idx, rxn_equation in enumerate(rxn_equations_list):
            if 'custom_energy_tuples' in kwargs:
                energy_tuple = \
                    self.get_relative_energy_tuple(kwargs['custom_energy_tuples'][idx])
            else:
                energy_tuple = self.get_energy_tuple(rxn_equation)
            x, y, ts_scale = self.get_potential_energy_points(
                energy_tuple=energy_tuple, n=n,
                subsection_length=subsection_length
            )
            x_scale = 2*subsection_length + ts_scale
            #get total y
            offseted_y = y + y_offset
            total_y.extend(offseted_y.tolist())
            #get total x
            offseted_x = x + x_offset
            total_x.extend(offseted_x.tolist())

            #######   collect values for adding annotations   #######

            #there are 3 points in each part
            p1 = (2*subsection_length+offseted_x[0], offseted_y[0])
            p3 = (offseted_x[-1], offseted_y[-1])
            if len(energy_tuple) == 3:
                p2 = (2.5*subsection_length+offseted_x[0], energy_tuple[1]+y_offset)
                piece_points.append((p1, p2, p3))
            elif len(energy_tuple) == 2:
                piece_points.append((p1, p3))
            ts_scales.append(ts_scale)

            ################   Collection End   ################

            #update offset value for next loop
            x_offset += x_scale
            y_offset = offseted_y[-1]

        #Add extral line

        #add head horizontal line
        head_extral_x = np.linspace(0.0, subsection_length, n).tolist()
        offseted_total_x = (np.array(total_x)+subsection_length).tolist()
        head_extral_y = [total_y[0]]*n
        total_x = head_extral_x + offseted_total_x
        total_y = head_extral_y + total_y
        #add the last horizontal line
        tail_extral_x = np.linspace(total_x[-1],
                                    total_x[-1]+subsection_length,
                                    n).tolist()
        tail_extral_y = [total_y[-1]]*n
        total_x += tail_extral_x
        total_y += tail_extral_y
        total_x, total_y = np.array(total_x), np.array(total_y)

        #get extreme values of x and y
        y_min, y_max = np.min(total_y), np.max(total_y)
        x_min, x_max = total_x[0], total_x[-1]
        y_scale = y_max-y_min
        ###################      start Artist Mode!      ######################

        multi = len(rxn_equations_list)/2.2
        fig = plt.figure(figsize=(multi*10, multi*5))
        #fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_ylim(y_min - y_scale/10, y_max + y_scale/10)
        ax.set_xlim(x_min-subsection_length, x_max+subsection_length)
        ax.set_title('Energy Profile',
                     fontdict={
                         'fontsize': 15,
                         'weight': 1000,
                         'verticalalignment': 'bottom'
                         })
        ax.set_xlabel('Reaction Coordinate')
        ax.set_ylabel('Free Energy(eV)')
        #add shadow
        if has_shadow:
            self.add_line_shadow(ax, total_x, total_y, depth=8,
                                 color='#595959')
        #main line
        ax.plot(total_x, total_y, color=line_color, linewidth=3)

        #get species list like [('COO^*(s)', ), ('2H^*(s)', 'H-H^*(s) +  *(s)')]
        tex_state_list = []
        for rxn_list in rxns_list:
            sub_tex_state_list = []
            for state_list in rxn_list[:-1]:
                tex_state = self.state_latex_convert(state_list)
                sub_tex_state_list.append(tex_state)
            tex_state_list.append(tuple(sub_tex_state_list))

        ####################     Main loop to add notes     ###################
        #This is the main loop to add extral info to main line
        #For each single loop, we get 3 or 2 points(tuple) for each elementary rxn
        #and the first 2 or 1 latex reaction state string(involved in list)

        def add_state_note(pts, tex_states):
            "function to add notes to a state"
            start_point, end_point = pts[0], pts[-1]

            #add horizontal lines and vertical arrowa
            if show_aux_line:
                #plot horizontal lines
                hori_x = np.linspace(start_point[0],
                                     end_point[0]+2*subsection_length, 50)
                hori_y = np.linspace(start_point[-1], start_point[-1], 50)
                ax.plot(hori_x, hori_y, color=line_color, linewidth=1,
                        linestyle='dashed')
                #plot vertical arrows
                if len(pts) == 3:
                    mid_point = pts[1]
                    #plot activation energy arrow
                    ax.annotate(
                        '',
                        xy=(mid_point[0], mid_point[-1]),
                        xycoords='data',
                        xytext=(mid_point[0], start_point[-1]),
                        textcoords='data',
                        arrowprops=dict(arrowstyle="<->")
                        )
                #plot reaction energy arrow
                ax.annotate(
                    '',
                    xy=(end_point[0]+subsection_length, end_point[-1]),
                    xycoords='data',
                    xytext=(end_point[0]+subsection_length, start_point[-1]),
                    textcoords='data',
                    arrowprops=dict(arrowstyle="<->")
                    )
            if show_note:
                #add species notes
                if len(pts) == 3:  # those who have TS
                    for state_idx, (pt, tex_state) in enumerate(zip(pts, tex_states)):
                        if state_idx == 0:  # for initial or final state
                            note_x = pt[0] - 1.8*subsection_length
                            note_y = pt[-1] + y_scale/50
                            ax.text(note_x, note_y, r'$\bf{'+tex_state+r'}$',
                                    fontdict={'fontsize': 13, 'color': '#1874CD'})
                        if state_idx == 1:  # for transition state
                            note_x = pt[0] - subsection_length/2
                            note_y = pt[-1] + y_scale/50
                            ax.text(note_x, note_y, r'$\bf{'+tex_state+r'}$',
                                    fontdict={'fontsize': 13, 'color': '#CD5555'})
                if len(pts) == 2:  # no TS
                    pt, tex_state = pts[0], tex_states[0]
                    note_x = pt[0] - 1.8*subsection_length
                    note_y = pt[-1] + y_scale/50
                    ax.text(note_x, note_y, r'$\bf{'+tex_state+r'}$',
                            fontdict={'fontsize': 13, 'color': '#1874CD'})

            if show_arrow:
                #add energy annotations
                el = Ellipse((2, -1), 0.5, 0.5)
                #annotation between IS and FS
                rxn_energy = round(end_point[-1] - start_point[-1], 2)
                rxn_energy_latex = r'$\bf{\Delta G = ' + str(rxn_energy) + r' eV}$'
                ax.annotate(rxn_energy_latex, xy=(end_point[0]+subsection_length,
                                                  (start_point[-1]+end_point[-1])/2),
                            xycoords='data',
                            xytext=(-70, 30), textcoords='offset points',
                            size=13, color='#8E388E',
                            arrowprops=dict(arrowstyle="simple",
                                            fc="0.6", ec="none",
                                            patchB=el,
                                            connectionstyle="arc3,rad=0.2"),
                            )
                if len(pts) == 3:
                    #annotation between TS and IS
                    act_energy = round((mid_point[-1] - start_point[-1]), 2)
                    act_energy_latex = r'$\bf{G_{a} = '+str(act_energy)+r' eV}$'
                    ax.annotate(act_energy_latex,
                                xy=(mid_point[0], (mid_point[1]+start_point[1])/2),
                                xycoords='data', xytext=(30, 20),
                                textcoords='offset points',
                                size=13, color='#FF4500',
                                arrowprops=dict(arrowstyle="simple",
                                                fc="0.6", ec="none",
                                                patchB=el,
                                                connectionstyle="arc3,rad=-0.2"),
                                )
        #########    Main Loop END    ########

        #use multiple thread to add notes
        note_threads = []
        nstates = len(tex_state_list)

        for pts, tex_states in zip(piece_points, tex_state_list):
            t = NoteThread(add_state_note, (pts, tex_states))
            note_threads.append(t)

        for i in xrange(nstates):
            note_threads[i].start()

        for i in xrange(nstates):
            note_threads[i].join()

        #add species note at last section
        if show_note:
            pt = pts[-1]
            tex_state = self.state_latex_convert(rxns_list[-1][-1])
            ax.text(pt[0]+0.2*subsection_length, pt[-1]+y_scale/50,
                    r'$\bf{'+tex_state+r'}$',
                    fontdict={'fontsize': 13, 'color': '#1874CD'})
        #remove xticks
        ax.set_xticks([])

    #####################   Artist Mode End   #####################

        #if kwargs.has_key('fname'):
        if 'fname' in kwargs:
            fname = kwargs['fname'] + '.' + fmt
        else:
            fname = 'multi_energy_diagram.' + fmt

        if show_mode == 'show':
            fig.show()
        elif show_mode == 'save':
            fig.savefig(fname)
        else:
            raise ValueError('Unrecognized show mode parameter : ' + show_mode)

        return fig, total_x, total_y
