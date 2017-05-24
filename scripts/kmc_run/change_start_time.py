#!/usr/bin/env python
# -*- coding: utf-8 -*-
globs, locs = {}, {}

exec(open("auto_lattice_trajectory_old.py", "r").read(), globs, locs)
times = locs["times"]

start_time = times[-1]

with open("pt-100.mkm", "r") as f:
    content = f.readlines()

found = False
for idx, line in enumerate(content):
    if line.startswith("start_time"):
        content[idx] = "start_time = {}".format(start_time)
        found = True

if not found:
    content.append('start_time = 0.0\n')

with open("pt-100.mkm", "w") as f:
    f.writelines(content)

