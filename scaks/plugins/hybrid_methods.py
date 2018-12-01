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
        model.logger().info('Use ODE integration to get new initial coverages...')

    end = 10
    span = 1e-2
    init_cvgs = model.coverages()

    new_cvgs = model.solver.solve_ode(time_end=end,
                                      time_span=span,
                                      initial_cvgs=init_cvgs)

    if model.log_allowed:
        model.logger().info('generate new initial coverages - success')

    return new_cvgs

