# KMC processes.
processes = [
    {
        "reaction": "CO_g + *_t -> CO_t",
        "description": "CO adsorption at top site.",
        "coordinates_group": [[[0.0, 0.0, 0.0], [-0.5, 0.0, 0.0],
                               [0.0, 0.5, 0.0], [0.5, 0.0, 0.0],
                               [0.0, -0.5, 0.0]]],
        "elements_before": ["V", "V", "V", "V", "V"],
        "elements_after": ["C", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "CO_g + *_b -> CO_b",
        "description": "CO adsorption at bridge site.",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, -0.5, 0.0]]],
        "elements_before": ["V", "V", "V"],
        "elements_after": ["C", "V", "V"],
        "basis_sites": [1],
    },
    {
        "reaction": "CO_g + *_b -> CO_b",
        "description": "CO adsorption at bridge site.",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]]],
        "elements_before": ["V", "V", "V"],
        "elements_after": ["C", "V", "V"],
        "basis_sites": [2],
    },
    {
        "reaction": "O2_g + 2*_b -> 2O_b",
        "description": "O2 dissociative adsorption at bridge sites basis site 1",
        "coordinates_group":[[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                             [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]],
        "elements_before": ["V", "V"],
        "elements_after": ["O_s", "O_s"],
        "basis_sites": [1],
    },
    {
        "reaction": "O2_g + 2*_b -> 2O_b",
        "description": "O2 dissociative adsorption at bridge sites basis site 2",
        "coordinates_group":[[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                             [[0.0, 0.0, 0.0], [0.0, -1.0, 0.0]]],
        "elements_before": ["V", "V"],
        "elements_after": ["O_s", "O_s"],
        "basis_sites": [2],
    },
    {
        "reaction": "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b",
        "description": "CO and O couple and desorption.",
        "coordinates_group":[[[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]],
                             [[0.0, 0.0, 0.0], [0.5, -0.5, 0.0]],
                             [[0.0, 0.0, 0.0], [-0.5, 0.5, 0.0]],
                             [[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0]]],
        "elements_before": ["V", "V"],
        "elements_after": ["O_s", "O_s"],
        "basis_sites": [1, 2],
    }
]

