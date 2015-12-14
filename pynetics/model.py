import os
import sys
import inspect
import re
import logging
import logging.config

from functions import *
from errors.error import *


class KineticModel(object):
    """
    Main class for kinetic models.
    """
    def __init__(self, **kwargs):

        #set physical constants
        self._kB = 8.617332478e-5  # Boltzmann constant from NIST
        self._h = 4.135667516e-15  # Planck constant from NIST

        #set kinetic attrs
        self._attr_type_dict = {
            'adsorbate_names': tuple,
            'transition_state_names': tuple,
            'gas_names': tuple,
            'descriptor_names': tuple,
            'surface_names': tuple,
            'species_definitions': dict,
            'elementary_rxns': tuple,
            'setup_file': str,
            'rxn_expression': list
        }

        # set logger
        self.set_logger()

        # parse in keyword args
        for key in kwargs:
            if key in self._attr_type_dict:
                try:
                    val = self._attr_type_dict[key](kwargs[key])
                except:
                    raise AttributeError('Argument \''+key+'\' '
                                         'is in wrong type.')
                setattr(self, key, val)
            else:  # if key is not essential attr for model
                self.logger.warning("redundant keyword - [ %s ]", str(key))

        #set elementary parse regex(compiled)
        self.regex_dict = {}

        states_regex = re.compile(r'([^\<\>]*)(?:\<?\-\>)' +
                                  r'(?:([^\<\>]*)(?:\<?\-\>))?([^\<\>]*)')
        self.regex_dict['IS_TS_FS'] = [states_regex, ['IS', 'TS', 'FS']]

        species_regex = re.compile(r'(\d*)([^\_\+\*\<\>]+)_(\d*)(\w+)')
        self.regex_dict['species'] = \
            [species_regex, ['stoichiometry', 'name', 'site_number', 'site']]

        site_regex = re.compile(r'(\d*)(?:\*\_)(\w+)')
        self.regex_dict['empty_site'] = \
            [site_regex, ['stoichiometry', 'site']]

#        self.regex_dict['species_separator'] = \
#        			[r'(?:\A|\s*\+\s*|\s+)', []]

        self.hasdata = False

        # load setup file
        if hasattr(self, 'setup_file'):
            self.logger.info('setup file [ %s ] is found', self.setup_file)
            model_name = self.setup_file.rsplit('.', 1)[0]
            setattr(self, 'model_name', model_name)
            self.load(self.setup_file)
            self.logger.info('kinetic modeling...success!\n')
        else:
            self.logger.warning('setup file not read...')

    def set_parser(self, parser_name):
        """
        Import parser and set the instance of it as attr of model
        """
        # hacked from CatMap (catmap/model.py)
        # https://docs.python.org/2/library/functions.html#__import__
        basepath = os.path.dirname(
            inspect.getfile(inspect.currentframe()))
        if basepath not in sys.path:
            sys.path.append(basepath)
        # from  loggers import logger
        _module = __import__('parsers', globals(), locals())
        parser_instance = getattr(_module, parser_name)(owner=self)
        setattr(self, 'parser', parser_instance)
        self.logger.info('parser is set.')

    def set_logger(self):
        """
        Get logging.logger instance as logger of kinetic model.
        """
        logger = logging.getLogger('model')
        if os.path.exists('./logging.conf'):
            logging.config.fileConfig('./logging.conf')
        else:
            logger.setLevel(logging.INFO)
            # create handlers
            std_hdlr = logging.FileHandler('out.log')
            std_hdlr.setLevel(logging.DEBUG)
            console_hdlr = logging.StreamHandler()
            console_hdlr.setLevel(logging.INFO)
            # create formatter and add it to the handlers
            formatter = logging.Formatter('%(name)s   %(levelname)-8s %(message)s')
            std_hdlr.setFormatter(formatter)
            console_hdlr.setFormatter(formatter)
            # add the handlers to logger
            logger.addHandler(std_hdlr)
            logger.addHandler(console_hdlr)

        self.logger = logger

    def load(self, setup_file):
        """
        Load 'setup_file' into kinetic model by exec setup file
        and assigning all local variables as attrs of model.
        For tools, create the instances of tool classes and
        assign them as the attrs of model.
        """
        self.logger.info('Loading Kinetic Model...\n')

        defaults = dict(
            data_file='data.pkl',
            decimal_precision=100,
            tools=['parser', 'plotter'],
            parser='RelativeEnergyParser',
            table_maker='CsvMaker',
            solver='SteadyStateSolver',
            corrector='ThermodynamicCorrector',
            plotter='EnergyProfilePlotter',
        )

        self.logger.info('read in parameters...')

        #exec setup file set local variables as attrs of model
        globs = {}
        locs = defaults
        execfile(setup_file, globs, locs)

        # customize model tools
        if 'tools' in locs:
            if 'parser' not in locs['tools']:
                raise ParameterError('[ parser ] must be in tools.')
            self.tools = locs['tools']
            del locs['tools']

            self.logger.info('tools = %s', str(self.tools))

        # assign parser ahead to provide essential attrs for other tools
        self.logger.info('instantiate %s', str(locs['parser']))
        self.set_parser(locs['parser'])

        # assign other tools
        for key in locs.keys():
            # ignore tools which will be loaded later
            if key in self.tools:
                continue
            # check type of variables
            if key in self._attr_type_dict:
                #chech attr type
                if type(locs[key]) != self._attr_type_dict[key]:
                    try:
                        locs[key] = self._attr_type_dict[key](locs[key])
                        setattr(self, key, locs[key])
                    except:
                        raise ValueError('\''+key+'\' is in wrong type. ' +
                                         str(self._attr_type_dict[key]) +
                                         ' object is expected.')
            setattr(self, key, locs[key])
            self.logger.info('%s = %s', key, str(locs[key]))

        #use parser parse essential attrs for other tools
        #parse elementary rxns
        self.logger.info('Parsing elementary rxns...')
        if self.rxn_expressions:
            self.parser.parse_elementary_rxns(self.rxn_expressions)

        # load tools of model
        self.logger.info('instantiate model tools...')
        for key in self.tools:
            # auto-import classes (HACKED from CatMap)
            if key == 'parser':  # ignore parser which is loaded before
                continue
            try:
                if locs[key]:
                    if not key.endswith('s'):
                        pyfile = key + 's'
                    else:
                        pyfile = key
                    basepath = os.path.dirname(
                        inspect.getfile(inspect.currentframe()))
                    if basepath not in sys.path:
                        sys.path.append(basepath)
                    sublocs = {}
                    _temp = \
                        __import__(pyfile, globals(), sublocs, [locs[key]])
                    tool_instance = getattr(_temp, locs[key])(owner=self)
                    setattr(self, key, tool_instance)
                    self.logger.info('%s = %s', key, locs[key])
                else:
                    setattr(self, key, None)
                    self.logger.warning('%s is set to None.')
            except ImportError:
                raise AttributeError(key.capitalize()+' '+locs[key] +
                                     ' could not be imported. ' +
                                     'Ensure that the class ' +
                                     'exists and is spelled properly.')
            #HACK END
