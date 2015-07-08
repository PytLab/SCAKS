'''
    Data for energy profile plotting.
'''

#elementary reaction equations
rxn_equations = [
    'FFA + *_s -> FFA_s',
    'H + *_s -> H_s',
    'FFA_s + H_s <-> FFA-H_s + *_s -> FFAH_s + *_s',
    'H + *_s -> H_s',
    'FFAH_s + H_s <-> FFAH-H_s + *_s -> FFAHH_s + *_s',
    'FFAHH_s <-> TS_s -> RO_s',
]

#relative energies
multi_energy_tuples = [
    # IS,  TS,  FS
    #--------------
    [
        (0.0, -1.09),
        (-1.09, -1.58),
        (-1.58, -0.57, -1.49),
        (-1.49, -1.98),
        (-1.98, -1.41, -1.49),
        (-1.49, -1.16, -1.91),
    ],

    [
        (0.0, -1.09),
        (-1.09, -1.58),
        (-1.58, -0.57, -1.49),
        (-1.49, -1.98),
        (-1.98, -1.41, -1.49),
        (-1.49, -1.37, -1.74),
    ]

]
