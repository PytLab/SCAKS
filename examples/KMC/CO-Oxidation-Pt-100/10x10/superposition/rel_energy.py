# Energy barriers.
#Ga = [0.0, 0.0, 0.0, 0.0, 0.0, 0.39]

# Reaction energies.
#dG = [-1.92, -2.09, -1.29, -2.33, -3.48, -0.46]

Ga, dG = [], []

# CO_g + *_t -> CO_t
Ga.append(0.0)
dG.append(-1.92)

# CO_g + *_b -> CO_b
Ga.append(0.0)
dG.append(-2.09)

# O2_g + 2*_t -> O2_2t
Ga.append(0.0)
dG.append(-1.13)

# O2_2t + 2*_b <-> O-O_2t + 2*_b -> 2O_b + 2*_t
Ga.append(0.85)
dG.append(-1.08)

# O2_2t + CO_b <-> OCO-O_2t + *_b -> O_b + CO2_g + 2*_t
Ga.append(0.78)
dG.append(-1.65)

# O2_g + 2*_b -> 2O_b
Ga.append(0.07)
dG.append(-2.39)

# CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b
Ga.append(0.39)
dG.append(-0.46)

# CO_b + O_b <-> OC-O_2b -> CO2_g + 2*_b
Ga.append(0.56)
dG.append(0.002)

# CO_b + *_t <-> CO_t + *_b -> CO_b + *_t
#Ga.append(0.17)
#dG.append(0.0)

# O_b + *_t <-> O_t + *_b -> O_b + *_t
Ga.append(1.15)
dG.append(0.0)

# O2_g + *_b -> O_b + O_b
Ga.append(0.91)
dG.append(-1.82)

# CO_b + O_b <-> CO-O_t + *_b -> CO2_g + *_b + *_t
Ga.append(1.15)
dG.append(-0.28)

