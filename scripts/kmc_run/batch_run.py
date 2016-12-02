#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import commands

import numpy as np

def dict2setup(d):
    setup = ""
    for key, value in d.iteritems():
        value = "'{}'".format(value) if type(value) is str else value
        line = "{} = {}\n".format(key, value)
        setup += line
    return setup

filename = "pt-100.mkm"

species_pressures = {"CO_g": np.arange(0.01, 0.1, 0.001)}

species_name = "CO_g"

if "__main__" == __name__:
    setup_file = "./template/{}".format(filename)
    if os.path.exists(setup_file):
        globs, locs = {}, {}
        execfile(setup_file, globs, locs)

    for p in species_pressures[species_name]:
        dest = str(p)
        print("Create job {}".format(p))

        os.mkdir(dest)
        commands.getstatusoutput("cp ./template/* {}".format(dest))
        locs['species_definitions'][species_name]['pressure'] = p

        # Write new setup file.
        with open("./{}/{}".format(dest, filename), "w") as f:
            content = dict2setup(locs)
            f.write(content)

        print("Ok.")

