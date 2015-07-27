'''
    Data for energy profile plotting.
'''

#elementary reaction equations
multi_rxn_equations = [
    [
        'CO2 + * -> CO2*',
        'CO2* <-> CO-O_* -> CO*',
        'CO* <-> C-O* -> C*',
    ],

    [
        'CO2 + * -> CO2*',
        'CO2* <-> CO-O_* -> CO*',
        'CO* + H* <-> CO-H* -> COH*',
        'COH* <-> CH-O* -> CH*',
    ],
]

#relative energies
multi_energy_tuples = [
    # IS,  TS,  FS
    #--------------
    [
        (0.0, -0.34),
        (-0.34, -0.08, -1.28),
        (-1.28, 1.49, -0.45),
    ],
    [
        (0.0, -0.34),
        (-0.34, -0.08, -1.28),
        (-1.28, -0.26, -0.66),
        (-0.66, 0.98, -0.89),
    ],
]

#line colors
colors = ['#A52A2A', '#000000']  # '#A52A2A', '#000000', '#36648B', '#FF7256', '#008B8B', '#7A378B'
