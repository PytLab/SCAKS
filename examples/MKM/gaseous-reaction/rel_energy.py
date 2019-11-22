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
Ga.append(0.0) #0.4
dG.append(1.82 - 2.56)

#5 'OCH3_n + OH_b + H_o <-> OCH3_n + OH-H_b + *_o -> OCH3_g + H2O_g + *_n + *_b',
Ga.append(3.45 - 1.82)
dG.append(0.41 - 1.82)

#6 'CH4_g + OH_b <-> CH3-H_g + OH_b -> CH3_i + H2O_g'
Ga.append(2.64 - 0.41)
dG.append(0.72 - 0.41)

#1 'OCH3_i <-> OCH2-H_i -> HCHO_g + H_i',
Ga.append(1.16)
dG.append(0.02)

#2 'OCH3_i + O2_g <-> OCH2-H_i + O2_g -> HCHO_g + OOH_i',
Ga.append(2.91)
dG.append(-1.10)

#3 'HCHO_g + CH3_i <-> CH3_i + H-CHO_g -> CHO_i + CH4_g',
Ga.append(1.56)
dG.append(-0.64)

#4 'CHO_i + O2_g <-> CO-H_i + O2_g -> CO_g + OOH_i',
Ga.append(0.78)
dG.append(-1.37)

#5 'CHO_i <-> CO-H_i -> CO_g + H_i',
Ga.append(0.69)
dG.append(-0.26)

#6 '2CH3_i <-> CH3-CH3_i -> C2H6_g'
Ga.append(0.0)
dG.append(-2.06)

#7 '2H_i -> H2_g',
Ga.append(0.0)
dG.append(-3.41)

#8 'H_i + CH4_g <-> CH3-H_g + H_i -> H2_g + CH3_i'
Ga.append(1.49)
dG.append(-0.18)

#9 'H_i + C2H6_g <-> C2H5-H_g + H_i -> H2_g + C2H5_i', #9
Ga.append(1.35)
dG.append(-0.45)

#10 'O2 + H_i <-> OO-H_i -> OOH_i', #10 
Ga.append(0.97)
dG.append(-1.12)

#11 'CH3_i + OOH_i <-> CH3-OOH_i -> CH3OOH_g', #11
Ga.append(0.0)
dG.append(-1.28)

#12 'CH3OOH_g <-> CH3O-OH_g -> CH3O_i + OH_i', #12
Ga.append(0.32)
dG.append(0.32)

#13 'OH_i + H_i -> H2O_g', #13
Ga.append(0.0)
dG.append(-4.04)

#14 'C2H6_g + CH3_i <-> C2H5-H_g + CH3_i -> C2H5_i + CH4_g', #14 
Ga.append(1.79)
dG.append(-0.27)

#15 'C2H5_i -> C2H4_g + H_i', #15
Ga.append(0.0)
dG.append(0.68)

