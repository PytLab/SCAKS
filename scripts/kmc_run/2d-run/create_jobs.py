#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import commands
import shutil

import numpy as np

def dict2setup(d):
    setup = ""
    for key, value in d.items():
        value = "'{}'".format(value) if type(value) is str else value
        line = "{} = {}\n".format(key, value)
        setup += line
    return setup

# Gas pressures.
pCOs = np.linspace(0.01, 20, 10)
pO2s = np.linspace(0.01, 30, 10)

filename = "pt-100.mkm"

if "__main__" == __name__:
    setup_file = "./template/{}".format(filename)
    if os.path.exists(setup_file):
        globs, locs = {}, {}
        exec(open(setup_file, "rb").read(), globs, locs)

    job_dirs = ""
    for pO2 in pO2s:
        out_dir = "pO2-{}".format(pO2)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        for pCO in pCOs:
            job_dir = "{}/pCO-{}".format(out_dir, pCO)
            print("INFO: create job {}".format(job_dir))
            job_dirs += "{}\n".format(job_dir)

            if os.path.exists(job_dir):
                shutil.rmtree(job_dir)

            os.mkdir(job_dir)
            commands.getstatusoutput("cp -r ./template/* {}".format(job_dir))
            locs['species_definitions']['CO_g']['pressure'] = pCO
            locs['species_definitions']['O2_g']['pressure'] = pO2

            # Write new setup file.
            with open("./{}/{}".format(job_dir, filename), "w") as f:
                content = dict2setup(locs)
                f.write(content)

    with open("job_dirs.txt", "w") as f:
        f.write(job_dirs)

