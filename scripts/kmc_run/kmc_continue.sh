#!/usr/bin/env bash

python="/data/apps/python/Python-2.7.5/bin/python"
change_start_time="/data/apps/python/Python-2.7.5/scripts/change_start_time.py"

create_new_job()
{
    local newdir="continue"
    mkdir ${newdir}
    for f in $(find ./ -type f)
    do
        cp ${f} ${newdir}
    done
    cd ${newdir}
    find ./ -name "auto_*" ! -name "auto_lattice_trajectory.py" | xargs rm
    mv auto_lattice_trajectory.py auto_lattice_trajectory_old.py
    echo "from auto_lattice_trajectory_old import types as t;types = t[-1]" > kmc_configuration.py
    rm *.o* *.script
    ${python} ${change_start_time}
}

main()
{
    create_new_job
}

main

