#!/usr/bin/env python
# -*- coding: utf-8 -*-

def ODE_integration(model, N):
    ''' Hybrid method using ODE integration

    :param model: current micro-kinetic model
    :type model: :obj:`scaks.models.MicroKineticModel`

    :param N: ODE integration attemp times
    :type N: int
    '''
    if model.log_allowed:
        model.logger.info('Use ODE integration to get new initial coverages...')

    end = 10**N
    span = 1e-2#10**-(N+2)#1e-2
    init_cvgs = model.solver.coverages
    model._MicroKinetcModel____ode_output_interval = int(end/span)//10

    new_cvgs = model.solver.solve_ode(time_end=end,
                                      time_span=span,
                                      initial_cvgs=init_cvgs)[-1]

    if model.log_allowed:
        model.logger.info('generate new initial coverages - success')

    return new_cvgs

