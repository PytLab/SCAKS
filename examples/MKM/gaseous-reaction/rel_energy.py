# Relative Energies for all elementary reactions.
Ga, dG = [], []

#1 'H2O_g + *_o + *_b <-> H-OH_b + *_o -> H_o + OH_b',
Ga.append(0.67)
dG.append(-0.25)

#2 'O2_g + *_n -> O2_n',
Ga.append(0.0)
dG.append(2.27)

#3 'O2_n + CH4_g + *_b <-> CH3-H_g + O2_n -> CH3_i + O_n + OH_b',
Ga.append(4.56 - 2.27)
dG.append(2.56 - 2.27)

#4 'CH3_i + O_n <-> CH3-O_n -> OCH3_n',
Ga.append(0.4)
dG.append(1.82 - 2.56)

