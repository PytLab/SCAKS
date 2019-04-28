# KMC processes.
processes = [
    # CO adsorbed at top site.
    {
        "reaction": "CH4_g + 2*_s <-> CH3-H_s + *_s ->  CH3_s + H_s",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],
             [0.5, 0.5, 0.0]],
        ]
        "elements_before": ['V', 'V', 'V'],
        "elements_after": ['V', 'C', 'H'],
        "basis_sites": [0],
    },

    {
        "reaction": "CH3_s + OH_s <-> CH3-OH_s + *_s -> CH3OH_g + 2*_s",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],  # CH3
             [0.5, 0.0, 0.0]], # OH
        ]

        "elements_before": ['V', 'C', 'O'],
        "elements_after": ['V', 'V', 'V'],
        "basis_sites": [0],
    },

]
