#elementary reaction equations
rxn_equations = [
    'HCOOH_g + 2*_s <-> HCOO-H_s + *_s -> HCOO_s + H_s',
    'HCOO_s + *_s <-> H-COO_s + *_s -> COO_s + H_s',
    'COO_s -> CO2_g + *_s',
    '2H_s <-> H-H_s + *_s -> H2_g + 2*_s'
]

#relative energies
energy_tuples = [
    #IS, TS, FS
    (0, 1, 0.5),
    (3, 4.7, 0.7),
    (0, 4),
    (3, 4.7, 0.7),
]
