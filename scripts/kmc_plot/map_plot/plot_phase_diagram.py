#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import interp2d
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt

from kynetix.compatutil import subprocess

# Get pO2 dirs.
pO2_dirs = [i for i in subprocess.getoutput('ls').split('\n') if i.startswith('pO2-')]

pO2s = []
tofs = []
cvgs_CO = []
cvgs_O = []

for pO2_dir in pO2_dirs:
    pO2 = float(pO2_dir.split('-')[-1])
    pO2s.append(pO2)
    # Get pCO dirs.
    cmd = "ls {}/".format(pO2_dir)
    pCO_dirs = [i for i in subprocess.getoutput(cmd).split('\n') if i.startswith('pCO-')]
    pCOs = [float(pCO_dir.split('-')[-1]) for pCO_dir in pCO_dirs]

    tofs_1d = []
    cvgs_CO_1d = []
    cvgs_O_1d = []

    for pCO_dir in pCO_dirs:
        print("Go to {}/{} ...".format(pO2_dir, pCO_dir))
        # Read tofs.
        filename = "{}/{}/auto_frequency.py".format(pO2_dir, pCO_dir)
        globs, locs = {}, {}
        exec(open(filename, "r").read(), globs, locs)
        reaction_rates = locs["reaction_rates"]
        tof = 0.0
        reactions = sorted(reaction_rates.keys())
        for idx in [0, 2, 8]:
            tof += reaction_rates[reactions[idx]]
        tofs_1d.append(tof)

        # Read coverages.
        filename = "{}/{}/auto_coverages.py".format(pO2_dir, pCO_dir)
        globs, locs = {}, {}
        exec(open(filename, "r").read(), globs, locs)
        cvg_O = locs["coverages"][-2][-1]
        cvg_CO = locs["coverages"][-1][-1]
        cvgs_O_1d.append(cvg_O)
        cvgs_CO_1d.append(cvg_CO)

    tofs.append(tofs_1d)
    cvgs_CO.append(cvgs_CO_1d)
    cvgs_O.append(cvgs_O_1d)

pCOs = np.array(pCOs)
pO2s = np.array(pO2s)
tofs = np.array(tofs)
cvgs_CO = np.array(cvgs_CO)
cvgs_O = np.array(cvgs_O)
cvgs = cvgs_O - cvgs_CO

# 2D interpolate.
pO2s.shape = (-1, 1)
interp_func_cvgs = interp2d(pCOs, pO2s, cvgs, kind="linear")
interp_func_tofs = interp2d(pCOs, pO2s, tofs, kind="linear")

# Plot 3D contour.
#y, x = np.mgrid[0:1:100j, 0:1:100j]
ynew = np.linspace(1e-5, 0.5, 100)
xnew = np.linspace(1e-5, 0.5, 200)
extent = [np.min(xnew), np.max(xnew), np.min(ynew), np.max(ynew)]
z_cvgs = interp_func_cvgs(xnew, ynew)
z_tofs = interp_func_tofs(xnew, ynew)

CS_cvgs = plt.contourf(xnew.reshape(-1), ynew.reshape(-1),
                       z_cvgs, 10, cmap=plt.cm.coolwarm)

CS_tofs = plt.contour(xnew.reshape(-1), ynew.reshape(-1), z_tofs, 10,
                      alpha=0.7, linewidths=1.0)

plt.xlabel("P(CO_g)/bar")
plt.ylabel(r"P(O_2_g)/bar")

# Make a colorbar.
cbar = plt.colorbar(CS_tofs)
cbar.ax.set_ylabel("TOF")

#cbar.add_lines(CS2)

plt.show()

