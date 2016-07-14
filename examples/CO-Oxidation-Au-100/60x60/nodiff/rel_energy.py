# Energy barriers.
#Ga = [0.0, 0.0, 0.0, 0.0, 0.0, 0.39]

# Reaction energies.
#dG = [-1.92, -2.09, -1.29, -2.33, -3.48, -0.46]

Ga, dG = [], []

# CO_g + *_t -> CO_t
Ga.append(0.0)
dG.append(-0.43)

# CO_g + *_b -> CO_b
Ga.append(0.0)
dG.append(-0.54)

# O2_g + 2*_b -> O2_2b
Ga.append(0.0)
dG.append(-1.02)

# O2_g + 2*_t -> O2_2t
Ga.append(0.0)
dG.append(-0.71)

# O2_2b <-> O-O_2b -> 2O_b
Ga.append(0.48)
dG.append(0.47)

# O2_2b + CO_b <-> O-OCO_3b -> O_b + CO2_g + 2*_b
Ga.append(0.42)
dG.append(-0.52)

# O_b + CO_t -> CO2_g + *_b + *_t
Ga.append(0.0)
dG.append(-3.14)

# CO_b + O_b -> CO2_g + 2*_b
Ga.append(0.0)
dG.append(-3.02)

# CO_b + O_b <-> O-CO_2b -> CO2_g + 2*_b
Ga.append(1.42)
dG.append(-3.01)

# CO_t + O_b <-> OC-O_b + *_t -> CO2_g + *_t + *_b
Ga.append(0.77)
dG.append(-2.05)

# CO_b + *_t <-> CO_t + *_b -> CO_b + *_t
#Ga.append(0.11)
#dG.append(0.0)

# O_b + *_t <-> O_t + *_b -> O_b + *_t
Ga.append(0.91)
dG.append(0.0)
