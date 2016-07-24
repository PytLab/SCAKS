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
dG.append(-0.56)

# O2_2t + 2*_f <-> 2*_f + O-O_2t -> 2O_f + 2*_t
Ga.append(0.98)
dG.append(-1.42)

# CO_t + O2_2t + *_f <-> OCO-O_3t + *_f -> CO2_g + O_f + 3*_t
Ga.append(0.45)
dG.append(-2.25)

# CO_g + O2_2t + *_t -> CO_t + O2_2t
Ga.append(0.0)
dG.append(-1.60)

# CO_t + O2_g + 2*_t -> CO_t + O2_2t
Ga.append(0.0)
dG.append(-0.59)

# CO_t + O2_2t + *_f <-> OCO-O_2t + *_t + *_f -> CO2_g + O_f + 3*_t
Ga.append(0.1)
dG.append(-2.39)

# O2_2t + CO_g + *_t -> CO_t + O2_2t
Ga.append(0.0)
dG.append(-1.46)

# O2_g + CO_t + 2*_t -> CO_t + O2_2t
Ga.append(0.0)
dG.append(-0.34)

# O_h + CO_f <-> O-CO_f + *_t -> CO2_g + *_t + *_f
Ga.append(0.51)
dG.append(-1.02)

# O_f + CO_h <-> O-CO_h + *_f -> CO2_g + *_t + *_f
Ga.append(0.82)
dG.append(-0.66)

# O_f + CO_t <-> O-CO_t + *_f -> *_t + *_f + CO2_g
Ga.append(0.87)
dG.append(-0.65)

# O_h + CO_t <-> O-CO_t + *_h -> *_t + *_h + CO2_g
Ga.append(0.43)
dG.append(-1.05)

# O2_g + *_h + *_f -> O_h + O_f
Ga.append(0.0)
dG.append(-1.05)

