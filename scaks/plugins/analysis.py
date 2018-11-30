#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .metaclasses import AnalysisMeta


class OnTheFlyAnalysis(metaclass=AnalysisMeta):
    ''' Class for providing an interface to easily extend and customize the behavior
    of the on-the-fly analysis in Newton iteration.

    Attribute:

        interval(:obj:`int`): The analysis interval in evolution iteration, default 
                              value is 1 meaning analyze every Newton step.
    '''
    # Analysis interval.
    interval = 1

    def setup(self, model, outer_counter):
        ''' Function called right before the start of a new Newton iteration to 
        allow for custom setup of the analysis object.

        :param model: current micro-kinetic model
        :type model: :obj:`scaks.models.MicroKineticModel`

        :param outer_counter: The outer count number (or the counter for current 
        initial coverages)
        :type outer_counter: int
        '''
        raise NotImplementedError

    def register_step(self, model, inner_counter, outer_counter):
        '''
        Function called in each iteration step.

        :param model: current micro-kinetic model
        :type model: :obj:`scaks.models.MicroKineticModel`

        :param inner_counter: The inner count number (or the counter for current 
        Newton step)
        :type inner_counter: int

        :param outer_counter: The outer count number (or the counter for current 
        initial coverages)
        :type outer_counter: int
        '''
        raise NotImplementedError

    def finalize(self, model, outer_counter):
        '''
        Called after each Newton iteration to allow for custom finalization and 
        post-processing of the collected data.

        :param model: current micro-kinetic model
        :type model: :obj:`scaks.models.MicroKineticModel`

        :param outer_counter: The outer count number (or the counter for current 
        initial coverages)
        :type outer_counter: int
        '''
        raise NotImplementedError

