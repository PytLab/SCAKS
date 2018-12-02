#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging


class AnalysisMeta(type):
    ''' Metaclass for analysis plugin class
    '''
    def __new__(cls, name, bases, attrs):
        # Check interval type.
        if 'interval' in attrs:
            interval = attrs['interval']
            if type(interval) is not int or interval <= 0:
                raise TypeError('analysis interval must be a positive integer')

        for method_name in ['setup', 'register_step', 'finalize']:
            method = attrs.get(method_name, None)
            if method is not None and not callable(method):
                msg = "{} must be a callable object".format(method)
                raise AttributeError(msg)
            # Set default interface methods.
            elif method is None:
                if method_name == 'setup':
                    attrs[method_name] = lambda self, model, outer_counter: None
                elif method_name == 'register_step':
                    attrs[method_name] = lambda self, model, inner_counter, outer_counter: None
                elif method_name == 'finalize':
                    attrs[method_name] = lambda self, model, outer_counter: None

        # Set logger.
        logger_name = 'scaks.{}'.format(name)
        attrs['logger'] = logging.getLogger(logger_name)

        return type.__new__(cls, name, bases, attrs)
