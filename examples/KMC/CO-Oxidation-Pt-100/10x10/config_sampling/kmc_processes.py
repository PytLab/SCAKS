# KMC processes.
processes = [
    # dissociation on adjacent top site.
    {
        "reaction": "O2_g + 2*_t -> 2O_t",
        "coordinates_group": [[[0.0, 0.0, 0.0],
                               [0.0, 0.5, 0.0],
                               [0.0, 1.0, 0.0]],

                              [[0.0, 0.0, 0.0],
                               [0.5, 0.0, 0.0],
                               [1.0, 0.0, 0.0]]
                             ],
        "elements_before": ['V', 'V', 'V'],
        "elements_after": ['O_s', 'V', 'O_s'],
        "basis_sites": [0]
    },
    # dissociation on adjacent bridge sites.
    {
        "reaction": "O2_g + 2*_b -> 2O_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]],
                              [[0.0, 0.0, 0.0], [0.5, -0.5, 0.0]],
                              [[0.0, 0.0, 0.0], [-0.5, 0.5, 0.0]],
                              [[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0]]],
        "elements_before": ['V', 'V'],
        "elements_after": ['O_s', 'O_s'],
        "basis_sites": [1],
    },

    # O2 dissociative adsorption.
    {
        "reaction": "O2_g + 2*_b -> O_b + O_b",
        "coordinates_group": [[[0.0, 0.0, 0.0],
                               [0.0, 0.5, 0.0],
                               [0.0, 1.0, 0.0]]],
        "elements_before": ['V', 'V', 'V'],
        "elements_after": ['O_s', 'V', 'O_s'],
        "basis_sites": [2],
    },
    {
        "reaction": "O2_g + 2*_b -> O_b + O_b",
        "coordinates_group": [[[0.0, 0.0, 0.0],
                               [0.5, 0.0, 0.0],
                               [1.0, 0.0, 0.0]]],
        "elements_before": ['V', 'V', 'V'],
        "elements_after": ['O_s', 'V', 'O_s'],
        "basis_sites": [1],
    },
    # O diffusion.
    {
        "reaction": "O_b + *_t <-> O_t + *_b -> O_b + *_t",
        "coordinates_group": [[[0.5, 0.0, 0.0],
                               [1.0, 0.0, 0.0],
                               [1.5, 0.0, 0.0]],

                              [[0.0, 0.5, 0.0],
                               [0.0, 1.0, 0.0],
                               [0.0, 1.5, 0.0]]],

        "elements_before": ['O_s', 'V', 'V'],
        "elements_after": ['V', 'V', 'O_s'],
        "basis_sites": [0],
    },
]

