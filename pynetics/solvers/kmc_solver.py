import logging
from math import exp

try:
    from KMCLib import *
except ImportError:
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "!!!                                                   !!!"
    print "!!!          WARNING: KMCLib is not installed         !!!"
    print "!!! Any kMC calculation using KMCLib will be disabled !!!"
    print "!!!                                                   !!!"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

from .. import KineticCoreComponent
from ..errors.error import *
from ..database.thermo_data import kB_eV
from ..database.lattice_data import *


class KMCSolver(KineticCoreComponent):
    def __init__(self, owner):
        '''
        Class for kinetic Monte Carlo simulation process.
        '''
        KineticCoreComponent.__init__(self, owner)

        # set logger
        self.logger = logging.getLogger('model.solvers.KMCSolver')

    def get_elementary_rate(self, elementary_rxn_list, force_TST=False):
        '''
        Function to get elementary reaction rate.

        Parameters:
        -----------
        elementary_rxn_list: elementary reaction states list, list of lists of str.

        force_TST: whether force to use TST to get rates, bool

        Returns:
        --------
        Rf: forward rate, float

        Rr: reversed rate, float.

        Example:
        --------
        >>> m.solver.get_elementary_rate([['CO_g', '*_s'], ['CO_s']])
        >>> (2262.3375403296886, 0.022775493982398507)

        '''
        self.logger.info('getting elementary reaction rates for %s',
                         str(elementary_rxn_list))
        # check input validity
        try:
            idx = self._owner.elementary_rxns_list.index(elementary_rxn_list)
        except ValueError:
            msg = '%s is not in elementary_rxns_list.' % str(elementary_rxn_list)
            raise ReactionEquationError(msg)

        # forced to use TST
        if force_TST:
            if not ('Ga' and 'dG' in self._owner.relative_energies):
                msg = 'No [ Ga ] and [ dG ] read, check your rel_energy.py please.'
                raise ParameterError(msg)
            # get energy info
            Ga = self._owner.relative_energies['Ga'][idx]
            dG = self._owner.relative_energies['dG'][idx]
            Gar = Ga - dG
            # get forward and reversed rates
            Rf = self.get_reaction_rate(Ga)
            Rr = self.get_reaction_rate(Gar)

        # use corresponding methods
        else:
            # get forward barrier
            Ea = self._owner.relative_energies['Ea'][idx]
            try:
                dE = self._owner.relative_energies['dG'][idx]
                free_energy = True
            except KeyError:
                dE = self._owner.relative_energies['dE'][idx]
                free_energy = False
            # forward rate
            Rf = self.get_forward_rate(elementary_rxn_list, Ea)
            # reversed rate (use dE as nonfree energy change)
            Rr = self.get_reversed_rate(elementary_rxn_list, Ea, dE,
                                        free_energy=free_energy)

        return Rf, Rr

    def check_gas_participating(self, reactants):
        '''
        Check whether gas species is in reactants.

        Parameters:
        -----------
        reactants: reactant species list, list of str.

        Returns:
        --------
        has_gas: gas participating or not, bool
        '''
        gas_name = self.extract_gas_name(reactants)

        if gas_name:
            return True
        else:
            return False

    def extract_gas_name(self, species_list):
        '''
        Extract gas name from a species list.
        '''
        # gas participating or not
        gases = []
        for sp in species_list:
            if sp.endswith('_g'):
                gas_name = self.split_species(sp)[-1]
                gases.append(gas_name)
        # check gas species number
        if len(gases) > 1:
            msg = ('There are more than one gases %s in species_list %s.' %
                   (str(gases), str(species_list)))
            raise ReactionEquationError(msg)
        # get gas name
        if gases:
            gas_name, = gases
        else:
            gas_name = None

        return gas_name

    def get_forward_rate(self, elementary_rxn_list, Ea):
        '''
        Function to determine the reaction type and get forward reaction rate.

        Parameters:
        -----------
        reactants: reactant species list, list of str.

        Ea: elementary reaction enengy barrier, float
        '''
        self.logger.info('getting forward rate for %s', str(elementary_rxn_list))
        reactants = elementary_rxn_list[0]

        gas_name = self.extract_gas_name(reactants)

        if gas_name:
            self.logger.info('%s is adsorption process, use Collision Theory.',
                             str(elementary_rxn_list))
            Rf = self.get_adsorption_rate(gas_name, Ea)
        else:
            self.logger.info('%s is not adsorption process, use TST.',
                             str(elementary_rxn_list))
            Rf = self.get_reaction_rate(Ea)

        return Rf

    def get_reversed_rate(self, elementary_rxn_list, Ea, dE, free_energy=False):
        '''
        Function to determine the reaction type and get reversed reaction rate.

        Parameters:
        -----------
        elementary_rxn_list: elementary reaction states list, list of lists of str.

        Ea: free energy barrier for forward reaction, float.

        dE: reaction energy change(free energy or not), float.

        free_energy: use dE as a free energy or not, bool.

        Returns:
        --------
        Rr: reversed rate of the elementary reaction, float.

        Example:
        --------
        >>> m.solver.get_reversed_rate([['CO_g', '*_s'], ['CO_s']], 0.0, -1.6)
        >>> 0.022775493982398507

        '''
        self.logger.info('getting reversed rate for %s', str(elementary_rxn_list))

        Ear = Ea - dE  # reversed reaction barrier
        reactants = elementary_rxn_list[0]

        is_desorption = self.check_gas_participating(reactants)
        if is_desorption:
            # use balance condition
            self.logger.info('%s is adsorption process, use balance condition.',
                             str(elementary_rxn_list))
            gas_name = self.extract_gas_name(reactants)
            Rr = self.get_desorption_rate(gas_name, dE, free_energy=free_energy)
        else:
            self.logger.info('%s is not adsorption process, reverse reaction equation ' +
                             'and get forward rate of it', str(elementary_rxn_list))
            reversed_rxn_list = list(reversed(elementary_rxn_list))
            Rr = self.get_forward_rate(reversed_rxn_list, Ear)

        return Rr

    def get_adsorption_rate(self, gas_name, Ea=0.0):
        '''
        Function to get gas adsorption rate using Collision Theory.

        Parameters:
        -----------
        gas_name: gas molecular formula with suffix, str.

        Ea: energy barrier (not free energy), float.

        Example:
        --------
        >>> m.solver.get_adsorption_rate('CO_g')
        >>> 2262.3375403296886

        '''
        # gas_name without suffix
        stripped_name = gas_name.split('_')[0]
        # parameters for collision theory
        act_ratio = self._owner.active_ratio                       # active area ratio
        Auc = self._owner.unitcell_area                            # unit cell area
        p = self._owner.species_definitions[gas_name]['pressure']  # partial pressure
        m = self._owner.parser.get_molecular_mass(stripped_name, absolute=True)  # molecular mass
        T = self._owner.temperature
        # forward rate
        Rads = self.get_kCT(Ea, Auc, act_ratio, p, m, T)

        return Rads  # s^-1

    def get_desorption_rate(self, gas_name, dE, Ea=0.0, free_energy=False):
        '''
        Function to get desorption in detailed balance condition.

        Parameters:
        -----------
        gas_name: gas molecular formula with suffix, str.

        Ea: energy barrier (not free energy), float.

        dE: correspinding adsorption energy change(E* - E_gas), float.

        free_energy: dE is free energy change or not, False by default, bool.

        Example:
        --------
        >>> m.solver.get_desorption_rate('CO_g', -1.6)
        >>> 0.022775493982398507
        '''
        Rads = self.get_adsorption_rate(gas_name, Ea)

        # gas_name without suffix
        stripped_name = gas_name.split('_')[0]

        T = self._owner.temperature

        if free_energy:
            self.logger.info('free energy read, use it directly.')
            Rdes = Rads/exp(-dE/(kB_eV*T))
        else:
            self.logger.info('NO free energy read, add correction automatically.')
            # get entropy contribution in gas free energy
            delta_miu = self._owner.corrector.entropy_correction(stripped_name)
            self.logger.info('%s entropy correction = %.3e', gas_name, delta_miu)
            K = exp((delta_miu - dE)/(kB_eV*T))  # equilibrum constant
            Rdes = Rads/K

        return Rdes

    def get_reaction_rate(self, Ga):
        '''
        Function to get reaction rate using Transition State Theory.

        Parameters:
        -----------
        Ea: energy barrier (not free energy), float.
        '''
        T = self._owner.temperature
        R = self.get_kTST(Ga, T)

        return R


class KMCLibSolver(KMCSolver):
    def __init__(self, owner):
        '''
        Class for kinetic Monte Carlo simulation process using KMCLib.
        '''
        KMCSolver.__init__(self, owner)

        # set logger
        self.logger = logging.getLogger('model.solvers.KMCLibSolver')

    def get_elementary_processes(self, elementary_rxn_list):
        '''
        Function to get KMCLib processes for an elementary reaction.
        '''
        rxn_expression = self.elementary_rxn_list2str(elementary_rxn_list)
        self.logger.info('getting process for [ %s ]', rxn_expression)

        # get center site and neighbors coordinates
        grid_type = self._owner.grid_type
        if grid_type not in grid_neighbor_offsets:
            raise GridTypeError('Unsupported grid type [ %s ]' % grid_type)
        neighbor_offsets = grid_neighbor_offsets[grid_type]
        coordinates = [(0.0, 0.0, 0.0)]
        coordinates.extend(neighbor_offsets)

        # get rate constants
        kf, kr = self.get_elementary_rate(elementary_rxn_list)

        # get elements changes
        get_elements_changes = self._owner.parser.get_elementary_elements_changes

        def get_single_direction_processes(elementary_rxn_list):
            ''' get single direction process objects '''
            elements_changes = get_elements_changes(elementary_rxn_list)
            processes = []
            for elements_change in elements_changes:
                elements_before, elements_after = elements_change
                self.logger.info('%s -> %s',
                                 str(elements_before),
                                 str(elements_after))
                p = KMCProcess(coordinates=coordinates,
                               elements_before=elements_before,
                               elements_after=elements_after,
                               basis_sites=[0],
                               rate_constant=kf)
                processes.append(p)

            return processes

        # forward direction
        self.logger.info('instantiating forward reaction processes...')
        fprocesses = get_single_direction_processes(elementary_rxn_list)

        # reversed direction
        self.logger.info('instantiating reversed reaction processes...')
        elementary_rxn_list.reverse()
        rprocesses = get_single_direction_processes(elementary_rxn_list)

        processes = fprocesses + rprocesses

        return processes
