#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from vaspy.atomco import PosCar

def get_site_coordinates(sites, origin, h, unit_bases):
    ''' 获取所有sites的在supercell中的分数坐标.

    Parameters:
    -----------
    sites: 所有吸附位点的相对坐标
    origin: 原点原子的分数坐标
    h: 所有site高度
    unit_bases: unitcell的基向量
    '''
    unit_bases = np.matrix(unit_bases).T
    sites = np.matrix(sites).T
    coords = unit_bases*sites
    origin = np.matrix(origin).T
    coords += origin
    coords[-1, :] += h

    # PDC.
    coords = np.array(coords)
    coords[coords > 1.0] -= 1.0

    return coords.T.tolist()

def get_config_poscar(config, positions, poscar_template='POSCAR'):
    poscar = PosCar(poscar_template)
    for atom, coord in zip(config, positions):
        if atom.startswith('O'):
            poscar.add_atom('O', coord)

    return poscar

if __name__ == '__main__':
    from auto_lattice_trajectory import sites
    origin = [0.201250, 0.158250, 0.370509]
    h = 0.08
    unit_bases = [[0.25, 0.0, 0.0],
                  [0.0, 0.25, 0.0],
                  [0.0, 0.0, 0.0]]
    site_coords = get_site_coordinates(sites, origin, h, unit_bases)

    from auto_lattice_trajectory import types
    for i, config in enumerate(types):
        poscar = get_config_poscar(config, site_coords)
        filename = 'configs/config_{}.cif'.format(i)
        with open(filename, 'w') as f:
            f.write(poscar.get_cif_content())
        print('{} created'.format(filename))

