# Energy barriers.
#Ga = [0.0, 0.0, 0.0, 0.0, 0.0, 0.39]

# Reaction energies.
#dG = [-1.92, -2.09, -1.29, -2.33, -3.48, -0.46]

Ga, dG = [], []

# CO_g + *_t -> CO_t
Ga.append(0.0)
dG.append(-1.58)

# CO_g + *_b -> CO_b
Ga.append(0.0)
dG.append(-1.68)

# CO_g + *_h -> CO_h
Ga.append(0.0)
dG.append(-1.71)

# CO_g + *_f -> CO_f
Ga.append(0.0)
dG.append(-1.75)

# O2_g + 2*_t -> O2_2t
Ga.append(0.0)
dG.append(0.56)

# O2_2t + 2*_f <-> 2*_f + O-O_2t -> 2O_f + 2*_t
Ga.append(0.98)
dG.append(-1.42)

