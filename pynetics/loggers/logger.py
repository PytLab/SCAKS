from string import Template
import cPickle


class Logger(object):
    """
    A general logger class for model and model's tools.
    """
    def __init__(self, owner):
        self._owner = owner
        self._templates_dict = {
            'attribute_load_failed': 'no attribute named ' +
                                     '${attribute} in setup file'
        }
        self._event_lines = []
        self._iter_lines = []
        self._warnings = []

        #attrs for archive data
        self.data_dict = {}  # object to be serialized
        self.data_file = 'data.pkl'  # pickle file name

    @staticmethod
    def write_logfile(filename, line):
        f = open(filename, 'a')
        f.write(line)
        f.close()

    def log(self, log_type, event, **kwargs):
        #event must in "operation_status" format
        message = self._templates_dict[event]
        operation, status = event.rsplit('_', 1)
        kwargs['operation'] = operation
        kwargs['status'] = status

        if status == 'warning':
            message = '${operation} warning! : ' + message
            message = Template(message).substitute(kwargs)
            #self._warnings.append(message)
            self.write_logfile('event.log', message)
        else:
            #iteration log
            if log_type == 'iteration':
                message = '${operation}_operation -> ${status} : ' + message
                message = Template(message).substitute(kwargs)
                #self._iter_lines.append(message)
                self.write_logfile('event.log', message+'\n')
            #log events
            if log_type == 'event':
                message = '${operation}_operation -> ${status} : ' + message
                message = Template(message).substitute(kwargs)
                #self._event_lines.append(message)
                self.write_logfile('event.log', message)
        print message

        return message

    def archive_data(self, data_name, data):
        "Update data dict and dump it to data file."
        #update data dict
        if data_name in self._owner.archived_variables:
            self.data_dict[data_name] = data
            #dump data dict to data file
            if self.data_dict:
                with open(self.data_file, 'wb') as f:
                    cPickle.dump(self.data_dict, f)

    def log_sscvg(self, cvgs_tuple, ads_names):
        "Log steady state coverage of every species."
        head_str = "\n %-5s     %-20s     %-30s\n" % \
                   ("Index", "Intermediate Name", "Steady State Coverage")
        line_str = '-'*60 + '\n'

        all_data = ''
        all_data += head_str + line_str
        for idx, (ads_name, cvg) in enumerate(zip(ads_names, cvgs_tuple)):
            idx = str(idx).zfill(2)
            data = " %-5s     %-20s     %-30.5e\n" % (idx, ads_name, float(cvg))
            all_data += data
        all_data += line_str

        print all_data
        self.write_logfile('event.log', all_data)

        return all_data

    def log_latex(self, latex_tup):
        "Append latex strings to 'formulas.tex'."
        latex_str = ''.join(latex_tup)
        latex_str += '\n'

        self.write_logfile('formulas.tex', latex_str)
        return latex_str
