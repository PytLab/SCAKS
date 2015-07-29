'''
    Data for energy profile plotting.
'''

#for a number of rxn equations
#COMMENT THEM OUT if not a number of rxn equations
'''
#elementary reaction equations
rxn_equations = [
    'HCOOH_g + 2*_s <-> HCOO-H_s + *_s -> HCOO_s + H_s',
    'HCOO_s + *_s <-> H-COO_s + *_s -> COO_s + H_s',
    'COO_s -> CO2_g + *_s',
    '2H_s <-> H-H_s + *_s -> H2_g + 2*_s'
]

#relative energies
energy_tuples = [
    # IS,  TS,  FS
    #--------------
    (0.0, 1.0, 0.5),
    (3.0, 4.7, 0.7),
    (0.0, 4.0),
    (3.0, 4.7, 0.7),
]
'''

#for a single rxn equation.
#COMMENT THEM OUT if not single rxn equation
'''
rxn_equation = 'HCOOH_g + 2*_s <-> HCOO-H_s + *_s -> HCOO_s + H_s'

#relative energies
energy_tuple = (0.0, 0.33, -0.26)
'''

#line color
color = '#000000'
