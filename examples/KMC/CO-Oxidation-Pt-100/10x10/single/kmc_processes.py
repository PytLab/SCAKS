# KMC processes.
processes = [
    # CO adsorbed at top site.
    {
        "reaction": "CH4_g + 2*_f <-> CH3-H_2f -> CH3_f + H_f",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],
             [0.5, 0.5, 0.0]],
        ],
        "elements_before": ['V', 'V', 'V'],
        "elements_after": ['V', 'C', 'H'],
        "basis_sites": [0],
    },

    {
        "reaction": "H2O2_g + *_f -> H2O2_f",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],  # CH3
             [0.5, 0.0, 0.0]], # OH
        ],

        "elements_before": ['V', 'C', 'V'],
        "elements_after": ['V', 'C', 'H1'],
        "basis_sites": [0],
    },
	{
        "reaction": "CH3_f + H2O2_f <-> CH3-OHOH_f + *_f -> CH3OH_f + OH_f",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],  # CH3
             [0.5, 0.0, 0.0]], # OH
        ],

        "elements_before": ['V', 'C', 'H1'],
        "elements_after": ['V', 'C1', 'OH'],
        "basis_sites": [0],
    },
	{
        "reaction": "CH3OH_f + OH_f + H_f <-> CH3OH-OHH_f + 2*_f -> CH3OH_f + H2O_f + *_f",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],  # CH3
             [0.5, 0.0, 0.0]], # OH
        ],

        "elements_before": ['V', 'C', 'H1'],
        "elements_after": ['V', 'C1', 'H2O'],
        "basis_sites": [0],
    },
	{
        "reaction": "CH3OH_f  -> CH3OH_g + *_f",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],  # CH3
             [0.5, 0.0, 0.0]], # OH
        ],

        "elements_before": ['V', 'C1', 'H2O'],
        "elements_after": ['V', 'V', 'H2O'],
        "basis_sites": [0],
    },
	{
        "reaction": "H2O_f -> H2O_g + *_f",
        "coordinates_group": [
            [[0.0, 0.0, 0.0],
             [0.0, 0.5, 0.0],  # CH3
             [0.5, 0.0, 0.0]], # OH
        ],

        "elements_before": ['V', 'V', 'H2O'],
        "elements_after": ['V', 'V', 'V'],
        "basis_sites": [0],
    },

]
