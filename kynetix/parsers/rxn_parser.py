# -*- coding:utf-8 -*-
import re

from kynetix.errors.error import *


class RxnEquation(object):
    """
    Class to create reaction equation object.
    """
    def __init__(self, rxn_equation):
        """
        Constructor.
        """
        self.__rxn_equation = rxn_equation

    def tolist(self):
        """
        Convert rxn_equation string to rxn_list(chem_state objects).
        """
        states_regex = re.compile(r'([^\<\>]*)(?:\<?\-\>)' +
                                  r'(?:([^\<\>]*)(?:\<?\-\>))?([^\<\>]*)')
        m = states_regex.search(self.__rxn_equation)
        state_list = []
        for idx in range(1, 4):
            if m.group(idx):
                chem_state_obj = ChemState(m.group(idx).strip())
                state_list.append(chem_state_obj)

        return state_list

    def to_formula_list(self):
        """
        Function to get list of formulas for the reaction equation.
        """
        state_list = []

        for state in self.tolist():
            formula_list = []
            for formula in state.tolist():
                formula_list.append(formula)
            state_list.append(formula_list)

        return state_list

    def check_conservation(self):
        """
        Function to check reaction equation conservation.

        Returns:
        --------
        Conservative or not, bool.
        """
        state_list = self.tolist()

        first = state_list[0]
        for state in state_list[1:]:
            first.conserve(state)

        # If passed, return True.
        return True

    def texen(self):
        """
        Get tex string of a reaction equation.
        """
        state_list = self.tolist()
        tex_list = [state_obj.texen() for state_obj in state_list]

        if len(tex_list) == 3:
            tex_str = (tex_list[0] + r' \leftrightarrow ' + tex_list[1] +
                       r' \rightarrow ' + tex_list[-1])
        elif len(tex_list) == 2:
            tex_str = tex_list[0] + r' \rightarrow ' + tex_list[-1]

        return tex_str

    # ----------------------------------------
    # Query functions.

    def rxn_equation(self):
        """
        Query function for reaction equation string.
        """
        return self.__rxn_equation


class ChemState(object):
    """
    Class to generate chemical state object.
    """
    def __init__(self, chem_state):
        """
        Constructor.
        """
        self.__chem_state = chem_state

    def split(self):
        """
        Function to split state to formula string list.
        """
        return [sp.strip() for sp in self.__chem_state.split('+')]

    def tolist(self):
        """
        Function to split state string to chemical formula list.
        """
        formula_str_list = self.split()
        formula_list = [ChemFormula(formula) for formula in formula_str_list]
        return formula_list

    def get_species_site_list(self):
        """
        Function to get species_site list of the state.
        """
        formula_list = self.tolist()

        species_site_list = [formula.species_site() for formula in formula_list]

        return species_site_list

    def get_species_site_dict(self):
        """
        Function to get species_site dictionary of the state.
        """
        formula_list = self.tolist()

        species_site_dict = {formula.species_site(): formula.stoichiometry()
                             for formula in formula_list}

        return species_site_dict

    def get_elements_dict(self):
        """
        Function to get element dictionary of the state.
        """
        # Get elements dict of all species.
        formula_list = self.tolist()
        elements_dicts = (formula.get_elements_dict() for formula in formula_list)

        # Merge all elements dicts.
        merged_dict = {}
        for elements_dict in elements_dicts:
            for element, number in elements_dict.iteritems():
                if element not in merged_dict:
                    merged_dict.setdefault(element, number)
                else:
                    merged_dict[element] += number

        return merged_dict

    def get_sites_dict(self):
        """
        Function to get sites dictionary of the state.
        """
        # Get sites dict of all formulas.
        formula_list = self.tolist()
        sites_dicts = (formula.get_sites_dict() for formula in formula_list)

        # Site types that are not included.
        exclusive_sites = ("g", "l")

        # Merge all elements dicts.
        merged_dict = {}
        for sites_dict in sites_dicts:
            for site, number in sites_dict.iteritems():
                # If gas or liquid, pass.
                if site in exclusive_sites:
                    continue

                if site not in merged_dict:
                    merged_dict.setdefault(site, number)
                else:
                    merged_dict[site] += number

        return merged_dict

    def conserve(self, another):
        """
        Function to check state conservation.

        Parameters:
        -----------
        another: Another ChemState instance.

        Returns:
        --------
        Conservative or not, bool.
        """
        # Check input parameter.
        if not isinstance(another, self.__class__):
            msg = "Parameter another must be instance of ChemState."
            return ParameterError(msg)

        # Check elements.
        dict1 = self.get_elements_dict()
        dict2 = another.get_elements_dict()

        if dict1 != dict2:
            msg_temp = "Mass of chemical states {} and {} are not conservative."
            msg = msg_temp.format(self.chem_state(), another.chem_state())
            raise ValueError(msg)

        # Check sites.
        dict1 = self.get_sites_dict()
        dict2 = another.get_sites_dict()

        if dict1 != dict2:
            msg_temp = "Site of chemical states {} and {} are not conservatvie."
            msg = msg_temp.format(self.chem_state(), another.chem_state())
            raise ValueError(msg)

        # If all tests passed, return True.
        return True

    def texen(self):
        """
        Get tex string.
        """
        formula_list = self.tolist()
        first_sp = formula_list[0]
        tex_str = first_sp.texen()
        for formula in formula_list[1:]:
            tex_str += r' + ' + formula.texen()

        return tex_str

    # --------------------------------------
    # Query functions.

    def chem_state(self):
        """
        Query function for state string.
        """
        return self.__chem_state


class ChemFormulaError(Exception):
    "Exception raised for errors in the chemical formula."
    pass


class ChemFormula(object):
    """
    Class to generate chemical formula object.
    """
    def __init__(self, formula):
        """
        Constructor.
        """
        self.__formula = formula
        self.__formula_regex = re.compile(r'(\d*)(([\w\*-]*)_(\d*)([a-z\*]+))')
        self.__sp_regex = re.compile(r'([a-zA-Z\*])(\d*)')

        self.__split()

    def __add__(self, formula_inst):
        """
        Overload + operation function.
        """
        chem_state = self.formula + ' + ' + formula_inst.formula
        return ChemState(chem_state)

    def __split(self):
        """
        Private helper function to split whole formual to
        stoichiometry, species name, site number, site name.
        """
        m = self.__formula_regex.search(self.__formula)
        if not m:
            msg = 'Unexpected chemical formula: {}'.format(self.__formula)
            raise ChemFormulaError(msg)
        else:
            self.__stoich = int(m.group(1)) if m.group(1) else 1
            self.__species_site = m.group(2)
            self.__species = m.group(3)
            self.__site = m.group(5)
            self.__nsite = int(m.group(4)) if m.group(4) else 1

    def type(self):
        """
        Function to get species type:  'gas' | 'liquid' | 'adsorbate'
        """
        if self.__site == "g":
            return "gas"
        elif self.__site == "l":
            return "liquid"
        else:
            if "*" in self.__species_site:
                return "site"
            else:
                return "adsorbate"

    def get_elements_dict(self):
        """
        Function to get elements dictionary of formula.
        """
        sp_elements_dict = self.get_species_elements_dict()

        # Get total elements dict.
        elements_dict = {elem: self.__stoich*num
                         for elem, num in sp_elements_dict.iteritems()}

        return elements_dict

    def get_species_elements_dict(self, species=None):
        """
        Function to split elements of species to element dict.

        Parameters:
        -----------
        species: The species name.

        Returns:
        --------
        The element dict of the species.
        """
        if not species:
            species = self.__species

        # If it's a site, no element returned.
        if species == "*":
            return {}

        element_list = self.__sp_regex.findall(species)

        element_dict = {}
        for element, number in element_list:
            number = int(number) if number else 1
            if element not in element_dict:
                element_dict.setdefault(element, number)
            else:
                element_dict[element] += number

        return element_dict

    def get_sites_dict(self):
        """
        Function to get site dictionary of formula.
        """
        single_dict = dict([(self.__site, self.__nsite)])
        sites_dict = {site: num*self.__stoich for site, num in single_dict.iteritems()}

        return sites_dict

    def conserve(self, another):
        """
        Function to check conservation.

        Parameters:
        -----------
        another: Another ChemFormula instance.

        Returns:
        --------
        Conservative or not, bool.
        """
        # Check parameter type.
        if not isinstance(another, self.__class__):
            msg = "Parameter another must be instance of ChemFormula."
            raise ParameterError(msg)

        # Check elements.
        dict1 = self.get_elements_dict()
        dict2 = another.get_elements_dict()

        if dict1 != dict2:
            msg_temp = "Mass of chemical formula {} and {} are not conservative."
            msg = msg_temp.format(self.formula(), another.formula())
            raise ValueError(msg)

        # Check sites.
        dict1 = self.get_sites_dict()
        dict2 = another.get_sites_dict()

        if dict1 != dict2:
            msg_temp = "Site of chemical formula {} and {} are not conservatvie."
            msg = msg_temp.format(self.formula(), another.formula())
            raise ValueError(msg)

        # If all tests passed, return True.
        return True


    def __sub_texen(self, sub_species):
        """
        Private helper function to get tex string of sub-species.
        """
        tex_str = r''
        splited_tuples = self.__sp_regex.findall(sub_species)

        for element, n in splited_tuples:
            if n:
                tex_str += element + r'_{' + n + r'}'
            else:
                tex_str += element

        return tex_str

    def texen(self):
        """
        Get tex string of species.
        """
        tex_str = r''
        tex_str += str(self.__stoich) if self.__stoich != 1 else ''

        # Get species tex string.
        if '-' in self.__species:
            sub_species_list = self.__species.split('-')
            tex_str += '-'.join([self.__sub_texen(sub_sp) for sub_sp in sub_species_list])
        else:
            tex_str += self.__sub_texen(self.__species)

        tex_str += r'(' + self.__site + r')'

        return tex_str

    # -----------------------------------------------
    # Query functions.

    def formula(self):
        """
        Query function for formula string.
        """
        return self.__formula

    def stoichiometry(self):
        """
        Query function for stoichiometry number.
        """
        return self.__stoich

    def species_site(self):
        """
        Query function for string of species on site.
        """
        return self.__species_site

    def species(self):
        """
        Query function for species name.
        """
        return self.__species

    def nsite(self):
        """
        Query function for sites number.
        """
        return self.__nsite

    def site(self):
        """
        Query function for sites type.
        """
        return self.__site
