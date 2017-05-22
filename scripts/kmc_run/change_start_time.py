#!/usr/bin/env python
# -*- coding: utf-8 -*-
from auto_lattice_trajectory_old import times

start_time = times[-1]

with open("pt-100.mkm", "r") as f:
    content = f.readlines()

for idx, line in enumerate(content):
    if line.startswith("start_time"):
        content[idx] = "start_time = {}".format(start_time)

with open("pt-100.mkm", "w") as f:
    f.writelines(content)

