import os
import sys
import inspect
import re

from functions import *


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

        self._tools = \
            ['parser', 'table_maker', 'solver', 'corrector', 'plotter']
            #['parser', 'table_maker', 'solver', 'thermodynamics', 'plotter']

        #set logger
        self.set_logger()
        #update model's own logger template
        self.logger_template_dict = {
            'keyword_set_warning': 'redundant keyword ->' +
                                   ' \'${keyword}=${value}\'',
            }
        self.logger._templates_dict.update(self.logger_template_dict)

        #parse in keyword args
        for key in kwargs:
            if key in self._attr_type_dict:
                try:
                    val = self._attr_type_dict[key](kwargs[key])
                except:
                    raise AttributeError('Argument \''+key+'\' '
                                         'is in wrong type.')
                setattr(self, key, val)
            else:  # if key is not essential attr for model
                self.logger.log(log_type='event', event='keyword_set_warning',
                                keyword=key, value=kwargs[key])

        #set elementary parse regex(compiled)
        self.regex_dict = {}

        states_regex = re.compile(r'([^\<\>]*)(?:\<?\-\>)' +
                                  r'(?:([^\<\>]*)(?:\<?\-\>))?([^\<\>]*)')
        self.regex_dict['IS_TS_FS'] = [states_regex, ['IS', 'TS', 'FS']]

        species_regex = re.compile(r'(\d*)([^\_\+\*\<\>]+)_(\w+)')
        self.regex_dict['species'] = \
            [species_regex, ['stoichiometry', 'name', 'site']]

        site_regex = re.compile(r'(\d*)(?:\*\_)(\w+)')
        self.regex_dict['empty_site'] = \
            [site_regex, ['stoichiometry', 'site']]

#        self.regex_dict['species_separator'] = \
#        			[r'(?:\A|\s*\+\s*|\s+)', []]

        self.hasdata = False

        #load setup file
        if hasattr(self, 'setup_file'):
            model_name = self.setup_file.rsplit('.', 1)[0]
            setattr(self, 'model_name', model_name)
            self.load(self.setup_file)
            #self.set_table_maker(self.table_maker_name)

    def set_parser(self, parser_name):
        """
        Import parser and set the instance of it as attr of model
        """
        #The 'BLACK MAGIC' is hacked from CatMap (catmap/model.py)
        #https://docs.python.org/2/library/functions.html#__import__
        basepath = os.path.dirname(
            inspect.getfile(inspect.currentframe()))
        if basepath not in sys.path:
            sys.path.append(basepath)
        #from loggers import logger
        _module = __import__('parsers', globals(), locals())
        parser_instance = getattr(_module, parser_name)(owner=self)
        setattr(self, 'parser', parser_instance)

    def set_logger(self):
        """
        import logger and get an instance of Logger class
        """
        #The 'BLACK MAGIC' is hacked from CatMap (catmap/model.py)
        #https://docs.python.org/2/library/functions.html#__import__
        basepath = os.path.dirname(
            inspect.getfile(inspect.currentframe()))
        if basepath not in sys.path:
            sys.path.append(basepath)
        #from loggers import logger
        _module = __import__('loggers.logger', globals(), locals())
        logger_instance = getattr(_module, 'Logger')(owner=self)
        setattr(self, 'logger', logger_instance)

    def load(self, setup_file):
        """
        Load 'setup_file' into kinetic model by exec setup file
        and assigning all local variables as attrs of model.
        For tools, create the instances of tool classes and
        assign them as the attrs of model.
        """
        defaults = dict(
            data_file='data.pkl',
            decimal_precision=100,
            parser='CsvParser',
            table_maker='CsvMaker',
            solver='SteadyStateSolver',
            #solver='QuasiEquilibriumSolver',
            corrector='ThermodynamicCorrector',
            plotter='ThermoPlotter'
        )

        #exec setup file set local variables as attrs of model
        globs = {}
        locs = defaults
        execfile(setup_file, globs, locs)

        #assign parser ahead to provide essential attrs for other tools
        self.set_parser(locs['parser'])

        #assign other tools
        for key in locs.keys():
            #ignore tools which will be loaded later
            if key in self._tools:
                #setattr(self, 'table_maker_name', locs[key])
                continue
            #check type of variables
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

        #use parser parse essential attrs for other tools
        #parse elementary rxns
        if self.rxn_expressions:
            self.parser.parse_elementary_rxns(self.rxn_expressions)

        #load tools of model
        for key in self._tools:
            #black magic to auto-import classes
            #HACKED from CatMap
            if key == 'parser':  # ignore parser loaded before
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
                else:
                    setattr(self, key, None)
            except ImportError:
                raise AttributeError(key.capitalize()+' '+locs[key] +
                                     ' could not be imported. ' +
                                     'Ensure that the class ' +
                                     'exists and is spelled properly.')
            #HACK END

    def make_logfile(self):
        #attributes log
        log_line_list = []
        #model's attributes
        model_head = '#'*10 + ' '*3 + 'Info of Model' + ' '*3 + '#'*10
        log_line_list.append(model_head)
        for attr in self.__dict__:
            if (not callable(self.__dict__[attr]) and
                    not attr.startswith('_')):
                log_line = attr + ' = ' + repr(self.__dict__[attr])
                log_line_list.append(log_line)

        #tools' log
        event_line_list = []
        for tool in self._tools:
            tool_obj = getattr(self, tool)
            #attr log
            tool_head = '#'*10 + ' '*3 + 'Info of ' + tool + ' '*3 + '#'*10
            log_line_list.append(tool_head)
            for attr in getattr(self, tool).__dict__:
                if (not callable(tool_obj.__dict__[attr]) and
                        not attr.startswith('_')):
                    log_line = attr + ' = ' + \
                        repr(tool_obj.__dict__[attr])
                    log_line_list.append(log_line)

            #event log
            for event_type in ['_warnings', '_event_lines', '_iter_lines']:
                event_line_list.extend(getattr(tool_obj.logger, event_type))
        event_content = '\n'.join(event_line_list)
        log_content = '\n\n'.join(log_line_list)
        #create log file
        f = open('info.log', 'w')
        f.write(log_content)
        f.close()

        f = open('event.log', 'w')
        f.write(event_content)
        f.close()
