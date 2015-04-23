from string import Template


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
