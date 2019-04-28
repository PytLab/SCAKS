# Relative Energies for all elementary reactions.
Ga, dG = [], []

# CH4_g + 2*_f <-> CH3-H_2f -> CH3_f + H_f
Ga.append(1.31)
dG.append(1.15)

# H2O2_g + *_f -> H2O2_f
Ga.append(0.0)
dG.append(0.19)

# CH3_f + H2O2_f <-> CH3-OHOH_f + *_f -> CH3OH_f + OH_f
Ga.append(0.69)
dG.append(-2.45)

# CH3OH_f + OH_f + H_f <-> CH3OH-OHH_f + 2*_f -> CH3OH_f + H2O_f + *_f
Ga.append(0.0)
dG.append(-1.3)

# CH3OH_f  -> CH3OH_g + *_f
Ga.append(0)
dG.append(-0.09)

# H2O_f -> H2O_g + *_f
Ga.append(0)
dG.append(0.13)

