# -*- coding:utf-8 -*-
import re

from kynetix.errors.error import *


class RxnEquation(object):
    """
    Class to create reaction equation object.
    """
    def __init__(self, rxn_equation):
        self.__rxn_equation = rxn_equation

    def tolist(self):
        "Convert rxn_equation string to rxn_list(chem_state objects)."
        states_regex = re.compile(r'([^\<\>]*)(?:\<?\-\>)' +
                                  r'(?:([^\<\>]*)(?:\<?\-\>))?([^\<\>]*)')
        m = states_regex.search(self.__rxn_equation)
        state_list = []
        for idx in range(1, 4):
            if m.group(idx):
                chem_state_obj = ChemState(m.group(idx).strip())
                state_list.append(chem_state_obj)

        return state_list

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
        self.__chem_state = chem_state
        self.__sp_list = [sp.strip() for sp in chem_state.split('+')]

    def texen(self):
        """
        Get tex string.
        """
        first_sp = ChemFormula(self.__sp_list[0])
        tex_str = first_sp.texen()
        for sp in self.__sp_list[1:]:
            tex_str += r' + ' + ChemFormula(sp).texen()

        return tex_str

    # --------------------------------------
    # Query functions.

    def chem_state(self):
        """
        Query function for state string.
        """
        return self.__chem_state

    def species_list(self):
        """
        Query function for species list in state.
        """
        return self.__sp_list

class ChemFormulaError(Exception):
    "Exception raised for errors in the chemical formula."
    pass


class ChemFormula(object):
    """
    Class to generate chemical formula object.
    """
    def __init__(self, formula):
        self.__formula = formula
        self.__formula_regex = re.compile(r'(\d*)([\w\*-]*)_(\d*)([a-z\*]+)')
        self.__sp_regex = re.compile(r'([a-zA-Z\*])(\d*)')

        self.__stoich, self.__species, self.__nsite, self.__site = self.__split()

    def __split(self):
        """
        Split whole formual to
        stoichiometry, species name, site number, site name.
        """
        m = self.__formula_regex.search(self.__formula)
        if not m:
            raise ChemFormulaError('Unexpected chemical formula: %s' %
                                   self.formula)
        else:
            stoich = int(m.group(1)) if m.group(1) else 1
            species = m.group(2)
            site = m.group(4)
            nsite = int(m.group(3)) if m.group(3) else 1
            return stoich, species, nsite, site

    def split_species(self, species=None):
        """
        Function to split elements of species.

        Parameters:
        -----------
        species: The species name.

        Returns:
        --------
        The element dict of the species.
        """
        if not species:
            species = self.__species

        element_list = self.__sp_regex.findall(species)

        element_dict = {}
        for element, number in element_list:
            number = int(number) if number else 1
            if element not in element_dict:
                element_dict.setdefault(element, number)
            else:
                element_dict[element] += number

        return element_dict

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

    def __add__(self, formula_inst):
        chem_state = self.formula + ' + ' + formula_inst.formula
        return ChemState(chem_state)

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
        dict1 = self.split_species()
        dict2 = another.split_species()

        if dict1 != dict2:
            msg_temp = "Mass of chemical formula {} and {} are not conservative."
            msg = msg_temp.format(self.formula(), another.formula())
            raise ValueError(msg)

        # Check sites type.
        msg_temp = "Site of chemical formula {} and {} are not conservatvie."
        if self.__site != another.site():
            msg = msg_temp.format(self.formula(), another.formula())
            raise ValueError(msg)

        # Check site number.
        if self.__nsite != another.nsite():
            msg = msg_temp.format(self.formula(), another.formula())
            raise ValueError(msg)

        # If all tests passed, return True.
        return True

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

