# KMC processes.
processes = [
    # CO adsorbed at top site.
    {
        "reaction": "CO_g + *_t -> CO_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 0.5, 0.0],  # 1, 2
                               [0.5, 0.5, 0.0], [0.5, 0.0, 0.0],  # 3, 4
                               [0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 5, 6
                               [-0.5, -0.5, 0.0], [-0.5, 0.0, 0.0],  # 7, 8
                               [-0.5, 0.5, 0.0]]],  # 9
        "elements_before": ["V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["C", "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },

    # CO adsorbed on bridge site, basis site 1.
    {
        "reaction": "CO_g + *_b -> CO_b",
        "coordinates_group": [[[0.0, 0.0, 0.0],  # 1
                               [0.0, 0.5, 0.0], [0.5, 0.0, 0.0],  # 2, 3
                               [0.0, -0.5, 0.0], [-0.5, 0.0, 0.0]]],  # 4, 5
        "elements_before": ["V", "V", "V", "V", "V"],
        "elements_after": ["C", "V", "V", "V", "V"],
        "basis_sites": [1],
    },

    # CO adsorbed on bridge site, basis site 2.
    {
        "reaction": "CO_g + *_b -> CO_b",
        "coordinates_group": [[[0.0, 0.0, 0.0],  # 1
                               [0.0, 0.5, 0.0], [0.5, 0.0, 0.0],  # 2, 3
                               [0.0, -0.5, 0.0], [-0.5, 0.0, 0.0]]],  # 4, 5
        "elements_before": ["V", "V", "V", "V", "V"],
        "elements_after": ["C", "V", "V", "V", "V"],
        "basis_sites": [2],
    },

    # O2 adsorbed lying.
    {
        "reaction": "O2_g + 2*_t -> O2_2t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],  # 1, 2
                               [-0.5, 0.0, 0.0], [-0.5, 0.5, 0.0],  # 3, 4
                               [0.0, 0.5, 0.0], [0.5, 0.5, 0.0],  # 5, 6
                               [1.0, 0.5, 0.0], [1.5, 0.5, 0.0],  # 7, 8
                               [1.5, 0.0, 0.0], [1.5, -0.5, 0.0],  # 9, 10
                               [1.0, -0.5, 0.0], [0.5, -0.5, 0.0],  # 11, 12
                               [0.0, -0.5, 0.0], [-0.5, -0.5, 0.0],  # 13, 14
                               [0.5, 0.0, 0.0]],  # 15
                              ],
        "elements_before": ["V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V"],
        "elements_after": ["O_r", "O_l", "V", "V", "V",
                           "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_g + 2*_t -> O2_2t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0],  # 1, 2
                               [0.0, -0.5, 0.0], [0.5, -0.5, 0.0],  # 3, 4
                               [0.5, 0.0, 0.0], [0.5, 0.5, 0.0],  # 5, 6
                               [0.5, 1.0, 0.0], [0.5, 1.5, 0.0],  # 7, 8
                               [0.0, 1.5, 0.0], [-0.5, 1.5, 0.0],  # 9, 10
                               [-0.5, 1.0, 0.0], [-0.5, 0.5, 0.0],  # 11, 12
                               [-0.5, 0.0, 0.0], [-0.5, -0.5, 0.0],  # 13, 14
                               [0.0, 0.5, 0.0]],  # 15
                              ],
        "elements_before": ["V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V"],
        "elements_after": ["O_u", "O_d", "V", "V", "V",
                           "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },

    # O2 dissociation directly.
    {
        "reaction": "O2_2t + 2*_b <-> O-O_2t + 2*_b -> 2O_b + 2*_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],  # 1, 2
                               [1.5, 0.0, 0.0], [-0.5, 0.0, 0.0],  # 3, 4
                               [-0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 5, 6
                               [0.5, 0.5, 0.0], [1.0, 0.5, 0.0],  # 7, 8
                               [1.5, 0.5, 0.0], [2.0, 0.0, 0.0],  # 9, 10
                               [1.5, -0.5, 0.0], [1.0, -0.5, 0.0],  # 11, 12
                               [0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 13, 14
                               [-0.5, -0.5, 0.0], [-1.0, 0.0, 0.0],  # 15, 16
                               [0.5, 0.0, 0.0]],  # 17
                              ],
        "elements_before": ["O_r", "O_l", "V", "V", "V",
                            "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "O_s", "O_s", "V",
                           "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_2t + 2*_b <-> O-O_2t + 2*_b -> 2O_b + 2*_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0],  # 1, 2
                               [0.0, 1.5, 0.0], [0.0, -0.5, 0.0],  # 3, 4
                               [0.5, -0.5, 0.0], [0.5, 0.0, 0.0],  # 5, 6
                               [0.5, 0.5, 0.0], [0.5, 1.0, 0.0],  # 7, 8
                               [0.5, 1.5, 0.0], [0.0, 2.0, 0.0],  # 9, 10
                               [-0.5, 1.5, 0.0], [-0.5, 1.0, 0.0],  # 11, 12
                               [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.0],  # 13, 14
                               [-0.5, -0.5, 0.0], [0.0, -1.0, 0.0],  # 15, 16
                               [0.0, 0.5, 0.0]],  # 17
                              ],
        "elements_before": ["O_u", "O_d", "V", "V", "V",
                            "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "O_s", "O_s", "V",
                           "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },

    # O2 adsorbed lying beside bridge CO.
    {
        "reaction": "O2_g + 2*_t -> O2_2t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-0.5, 0.0, 0.0],  # 1, 2, 3
                               [-0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 4, 5
                               [0.5, 0.5, 0.0], [1.0, 0.5, 0.0],  # 6, 7
                               [1.5, 0.5, 0.0], [1.5, 0.0, 0.0],  # 8, 9
                               [1.5, -0.5, 0.0], [1.0, -0.5, 0.0],  # 10, 11
                               [0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 12, 13
                               [-0.5, -0.5, 0.0], [-1.0, 0.0, 0.0],  # 14, 15
                               [0.5, 0.0, 0.0]],  # 16
                              ],
        "elements_before": ["V", "V", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_r", "O_l", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_g + 2*_t -> O2_2t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0],  # 1, 2, 3
                               [1.5, 0.5, 0.0], [1.0, 0.5, 0.0],  # 4, 5
                               [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 6, 7
                               [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.0],  # 8, 9
                               [-0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 10, 11
                               [0.5, -0.5, 0.0], [1.0, -0.5, 0.0],  # 12, 13
                               [1.5, -0.5, 0.0], [2.0, -0.5, 0.0],  # 14, 15
                               [0.5, 0.0, 0.0]],  # 16
                              ],
        "elements_before": ["V", "V", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_r", "O_l", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_g + 2*_t -> O2_2t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.5, 0.0],  # 1, 2, 3
                               [0.5, 0.5, 0.0], [0.5, 0.0, 0.0],  # 4, 5
                               [0.5, -0.5, 0.0], [0.5, -1.0, 0.0],  # 6, 7
                               [0.5, -1.5, 0.0], [0.0, -1.5, 0.0],  # 8, 9
                               [-0.5, -1.5, 0.0], [-0.5, -1.0, 0.0],  # 10, 11
                               [-0.5, -0.5, 0.0], [-0.5, 0.0, 0.0],  # 12, 13
                               [-0.5, 0.5, 0.0], [0.0, 1.0, 0.0],  # 14, 15
                               [0.0, -0.5, 0.0]],  # 16
                              ],
        "elements_before": ["V", "V", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_d", "O_u", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_g + 2*_t -> O2_2t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, -1.5, 0.0],  # 1, 2, 3
                               [0.5, -1.5, 0.0], [0.5, -1.0, 0.0],  # 4, 5
                               [0.5, -0.5, 0.0], [0.5, 0.0, 0.0],  # 6, 7
                               [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 8, 9
                               [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.0],  # 10, 11
                               [-0.5, -0.5, 0.0], [-0.5, -1.0, 0.0],  # 12, 13
                               [-0.5, -1.5, 0.0], [0.0, -2.0, 0.0],  # 14, 15
                               [0.0, -0.5, 0.0]],  # 16
                              ],
        "elements_before": ["V", "V", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_d", "O_u", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },

    # CO adsorbed lying beside bridge O2.
    {
        "reaction": "CO_g + *_b -> CO_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0],  # 1, 2, 3
                               [1.5, 0.5, 0.0], [1.0, 0.5, 0.0],  # 4, 5
                               [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 6, 7
                               [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.0],  # 8, 9
                               [-0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 10, 11
                               [0.5, -0.5, 0.0], [1.0, -0.5, 0.0],  # 12, 13
                               [1.5, -0.5, 0.0], [2.0, -0.5, 0.0],  # 14, 15
                               [0.5, 0.0, 0.0]],  # 16
                              ],
        "elements_before": ["O_r", "O_l", "V", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_r", "O_l", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "CO_g + *_b -> CO_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-0.5, 0.0, 0.0],  # 1, 2, 3
                               [-0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 4, 5
                               [0.5, 0.5, 0.0], [1.0, 0.5, 0.0],  # 6, 7
                               [1.5, 0.5, 0.0], [1.5, 0.0, 0.0],  # 8, 9
                               [1.5, -0.5, 0.0], [1.0, -0.5, 0.0],  # 10, 11
                               [0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 12, 13
                               [-0.5, -0.5, 0.0], [-1.0, 0.0, 0.0],  # 14, 15
                               [0.5, 0.0, 0.0]],  # 16
                              ],
        "elements_before": ["O_r", "O_l", "V", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_r", "O_l", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "CO_g + *_b -> CO_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.5, 0.0],  # 1, 2, 3
                               [0.5, 0.5, 0.0], [0.5, 0.0, 0.0],  # 4, 5
                               [0.5, -0.5, 0.0], [0.5, -1.0, 0.0],  # 6, 7
                               [0.5, -1.5, 0.0], [0.0, -1.5, 0.0],  # 8, 9
                               [-0.5, -1.5, 0.0], [-0.5, -1.0, 0.0],  # 10, 11
                               [-0.5, -0.5, 0.0], [-0.5, 0.0, 0.0],  # 12, 13
                               [-0.5, 0.5, 0.0], [0.0, 1.0, 0.0],  # 14, 15
                               [0.0, -0.5, 0.0]],  # 16
                              ],
        "elements_before": ["O_d", "O_u", "V", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_d", "O_u", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },
    {
        "reaction": "CO_g + *_b -> CO_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, -1.5, 0.0],  # 1, 2, 3
                               [0.5, -1.5, 0.0], [0.5, -1.0, 0.0],  # 4, 5
                               [0.5, -0.5, 0.0], [0.5, 0.0, 0.0],  # 6, 7
                               [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 8, 9
                               [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.0],  # 10, 11
                               [-0.5, -0.5, 0.0], [-0.5, -1.0, 0.0],  # 12, 13
                               [-0.5, -1.5, 0.0], [0.0, -2.0, 0.0],  # 14, 15
                               [0.0, -0.5, 0.0]],  # 16
                              ],
        "elements_before": ["O_d", "O_u", "V", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_d", "O_u", "C", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [0],
    },

    # O2 dissociation with CO.
    {
        "reaction": "O2_2t + CO_b <-> OCO-O_2t + *_b -> O_b + CO2_g + 2*_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0],  # 1, 2, 3
                               [1.5, 0.5, 0.0], [1.0, 0.5, 0.0],  # 4, 5
                               [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 6, 7
                               [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.0],  # 8, 9
                               [-0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 10, 11
                               [0.5, -0.5, 0.0], [1.0, -0.5, 0.0],  # 12, 13
                               [1.5, -0.5, 0.0], [2.0, -0.5, 0.0],  # 14, 15
                               [0.5, 0.0, 0.0]],  # 16
                              ],
        "elements_before": ["O_r", "O_l", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "O_s"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_2t + CO_b <-> OCO-O_2t + *_b -> O_b + CO2_g + 2*_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-0.5, 0.0, 0.0],  # 1, 2, 3
                               [-0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 4, 5
                               [0.5, 0.5, 0.0], [1.0, 0.5, 0.0],  # 6, 7
                               [1.5, 0.5, 0.0], [1.5, 0.0, 0.0],  # 8, 9
                               [1.5, -0.5, 0.0], [1.0, -0.5, 0.0],  # 10, 11
                               [0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 12, 13
                               [-0.5, -0.5, 0.0], [-1.0, 0.0, 0.0],  # 14, 15
                               [0.5, 0.0, 0.0]],  # 16
                              ],
        "elements_before": ["O_r", "O_l", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "O_s"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_2t + CO_b <-> OCO-O_2t + *_b -> O_b + CO2_g + 2*_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.5, 0.0],  # 1, 2, 3
                               [0.5, 0.5, 0.0], [0.5, 0.0, 0.0],  # 4, 5
                               [0.5, -0.5, 0.0], [0.5, -1.0, 0.0],  # 6, 7
                               [0.5, -1.5, 0.0], [0.0, -1.5, 0.0],  # 8, 9
                               [-0.5, -1.5, 0.0], [-0.5, -1.0, 0.0],  # 10, 11
                               [-0.5, -0.5, 0.0], [-0.5, 0.0, 0.0],  # 12, 13
                               [-0.5, 0.5, 0.0], [0.0, 1.0, 0.0],  # 14, 15
                               [0.0, -0.5, 0.0]],  # 16
                              ],
        "elements_before": ["O_d", "O_u", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "O_s"],
        "basis_sites": [0],
    },
    {
        "reaction": "O2_2t + CO_b <-> OCO-O_2t + *_b -> O_b + CO2_g + 2*_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, -1.5, 0.0],  # 1, 2, 3
                               [0.5, -1.5, 0.0], [0.5, -1.0, 0.0],  # 4, 5
                               [0.5, -0.5, 0.0], [0.5, 0.0, 0.0],  # 6, 7
                               [0.5, 0.5, 0.0], [0.0, 0.5, 0.0],  # 8, 9
                               [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.0],  # 10, 11
                               [-0.5, -0.5, 0.0], [-0.5, -1.0, 0.0],  # 12, 13
                               [-0.5, -1.5, 0.0], [0.0, -2.0, 0.0],  # 14, 15
                               [0.0, -0.5, 0.0]],  # 16
                              ],
        "elements_before": ["O_d", "O_u", "C", "V", "V", "V", "V", "V",
                            "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V",
                           "V", "V", "V", "V", "V", "V", "V", "O_s"],
        "basis_sites": [0],
    },

    # O2 dissociative adsorption basis site 1.
    {
        "reaction": "O2_g + 2*_b -> 2O_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                               [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0],
                               [1.0, 0.5, 0.0], [1.5, 0.0, 0.0],
                               [1.0, -0.5, 0.0], [0.0, -0.5, 0.0],
                               [0.5, 0.0, 0.0]],
                              ],
        "elements_before": ["V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_s", "O_s", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },

    # O2 dissociative adsorption basis site 2.
    {
        "reaction": "O2_g + 2*_b -> 2O_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                               [0.0, -0.5, 0.0], [-0.5, 0.0, 0.0],
                               [-0.5, 1.0, 0.0], [0.0, 1.5, 0.0],
                               [0.5, 1.0, 0.0], [0.5, 0.0, 0.0],
                               [0.0, 0.5, 0.0]],
                              ],
        "elements_before": ["V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["O_s", "O_s", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [2],
    },

    # CO2_g associative desorption.
    {
        "reaction": "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.5, -0.5, 0.0],  # 1, 2
                               [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0],  # 3, 4
                               [0.0, -0.5, 0.0], [0.5, -1.0, 0.0],  # 5, 6
                               [1.0, -0.5, 0.0], [0.0, 0.5, 0.0]],  # 7, 8
                              ],
        "elements_before": ["C", "O_s", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },
    {
        "reaction": "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.5, -0.5, 0.0],  # 1, 2
                               [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0],  # 3, 4
                               [0.0, -0.5, 0.0], [0.5, -1.0, 0.0],  # 5, 6
                               [1.0, -0.5, 0.0], [0.0, 0.5, 0.0]],  # 7, 8
                              ],
        "elements_before": ["O_s", "C", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },
    {
        "reaction": "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0],  # 1, 2
                               [-0.5, 0.0, 0.0], [0.5, 0.0, 0.0],  # 3, 4
                               [0.0, -0.5, 0.0], [-0.5, -1.0, 0.0],  # 5, 6
                               [-1.0, -0.5, 0.0], [0.0, 0.5, 0.0]],  # 7, 8
                              ],
        "elements_before": ["O_s", "C", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },
    {
        "reaction": "CO_b + O_b <-> CO-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [-0.5, -0.5, 0.0],  # 1, 2
                               [-0.5, 0.0, 0.0], [0.5, 0.0, 0.0],  # 3, 4
                               [0.0, -0.5, 0.0], [-0.5, -1.0, 0.0],  # 5, 6
                               [-1.0, -0.5, 0.0], [0.0, 0.5, 0.0]],  # 7, 8
                              ],
        "elements_before": ["C", "O_s", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },

    # CO2_g associative desorption basis site 1.
    {
        "reaction": "CO_b + O_b <-> OC-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                               [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0],
                               [1.0, 0.5, 0.0], [1.5, 0.0, 0.0],
                               [1.0, -0.5, 0.0], [0.0, -0.5, 0.0],
                               [0.5, 0.0, 0.0]],
                              ],
        "elements_before": ["C", "O_s", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },
    {
        "reaction": "CO_b + O_b <-> OC-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                               [-0.5, 0.0, 0.0], [0.0, 0.5, 0.0],
                               [1.0, 0.5, 0.0], [1.5, 0.0, 0.0],
                               [1.0, -0.5, 0.0], [0.0, -0.5, 0.0],
                               [0.5, 0.0, 0.0]],
                              ],
        "elements_before": ["O_s", "C", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },

    # CO2_g associative desorption basis site 2.
    {
        "reaction": "CO_b + O_b <-> OC-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                               [0.0, -0.5, 0.0], [-0.5, 0.0, 0.0],
                               [-0.5, 1.0, 0.0], [0.0, 1.5, 0.0],
                               [0.5, 1.0, 0.0], [0.5, 0.0, 0.0],
                               [0.0, 0.5, 0.0]],
                              ],
        "elements_before": ["C", "O_s", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [2],
    },
    {
        "reaction": "CO_b + O_b <-> OC-O_2b -> CO2_g + 2*_b",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                               [0.0, -0.5, 0.0], [-0.5, 0.0, 0.0],
                               [-0.5, 1.0, 0.0], [0.0, 1.5, 0.0],
                               [0.5, 1.0, 0.0], [0.5, 0.0, 0.0],
                               [0.0, 0.5, 0.0]],
                              ],
        "elements_before": ["O_s", "C", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [2],
    },

#    # CO diffusion at basis site 1.
#    {
#        "reaction": "CO_b + *_t <-> CO_t + *_b -> CO_b + *_t",
#        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0],  # 1, 2
#                               [0.0, 1.5, 0.0], [0.5, 1.0, 0.0],  # 3, 4
#                               [0.5, -0.5, 0.0], [0.0, -0.5, 0.0],  # 5, 6
#                               [-0.5, 0.0, 0.0], [-0.5, 1.0, 0.0],  # 7, 8
#                               [0.0, 0.5, 0.0], [-0.5, 0.5, 0.0],  # 9, 10
#                               [0.5, 0.5, 0.0]],  # 11
#                              ],
#        "elements_before": ["C", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
#        "elements_after": ["V", "C", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
#        "basis_sites": [1],
#    },
#
#    # CO diffusion at basis site 2.
#    {
#        "reaction": "CO_b + *_t <-> CO_t + *_b -> CO_b + *_t",
#        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],  # 1, 2
#                               [1.5, 0.0, 0.0], [1.0, -0.5, 0.0],  # 3, 4
#                               [0.0, -0.5, 0.0], [-0.5, 0.0, 0.0],  # 5, 6
#                               [0.0, 0.5, 0.0], [1.0, 0.5, 0.0],  # 7, 8
#                               [0.5, 0.0, 0.0], [0.5, 0.5, 0.0],  # 9, 10
#                               [0.5, -0.5, 0.0]],  # 11
#                              ],
#        "elements_before": ["C", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
#        "elements_after": ["V", "C", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
#        "basis_sites": [2],
#    },

    # O diffusion at basis site 1.
    {
        "reaction": "O_b + *_t <-> O_t + *_b -> O_b + *_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0],  # 1, 2
                               [0.0, 1.5, 0.0], [0.5, 1.0, 0.0],  # 3, 4
                               [0.5, 0.0, 0.0], [0.0, -0.5, 0.0],  # 5, 6
                               [-0.5, 0.0, 0.0], [-0.5, 1.0, 0.0],  # 7, 8
                               [0.0, 0.5, 0.0], [-0.5, 0.5, 0.0],  # 9, 10
                               [0.5, 0.5, 0.0]],  # 11
                              ],
        "elements_before": ["O_s", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "O_s", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [1],
    },

    # O diffusion at basis site 2.
    {
        "reaction": "O_b + *_t <-> O_t + *_b -> O_b + *_t",
        "coordinates_group": [[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0],  # 1, 2
                               [1.5, 0.0, 0.0], [1.0, -0.5, 0.0],  # 3, 4
                               [0.0, -0.5, 0.0], [-0.5, 0.0, 0.0],  # 5, 6
                               [0.0, 0.5, 0.0], [1.0, 0.5, 0.0],  # 7, 8
                               [0.5, 0.0, 0.0], [0.5, 0.5, 0.0],  # 9, 10
                               [0.5, -0.5, 0.0]],  # 11
                              ],
        "elements_before": ["O_s", "V", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "elements_after": ["V", "O_s", "V", "V", "V", "V", "V", "V", "V", "V", "V"],
        "basis_sites": [2],
    },
]