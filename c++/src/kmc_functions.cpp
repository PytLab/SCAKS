#ifndef PLUGIN_BACKENDS_
#include "kmc_functions.h"
#endif

std::vector<double> collect_coverages(const std::vector<std::string> & types,
                                      const std::vector<std::string> & possible_types)
{
    // Get ntotal umber of sites.
    const size_t nsite = types.size();

    // Number of different types.
    std::vector<double> ntypes(possible_types.size(), 0.0);

    // Loop over all types to collect coverages.
    for (size_t type_idx = 0; type_idx < types.size(); ++type_idx)
    {
        // Current element.
        const std::string & element = types[type_idx];

        // Find the element in possible types, accumulate coverage.
        DoubleIterType ntypes_it = ntypes.begin();
        ConstStrIterType possible_types_it = possible_types.begin();
        const ConstStrIterType end = possible_types.end();

        for (; possible_types_it != end; ++ntypes_it, ++possible_types_it)
        {
            if (element == *possible_types_it)
            {
                // Accumulate.
                *ntypes_it += 1.0;
            }
        }
    }

    // Calculate coverages.
    std::vector<double> coverages;
    double coverage;

    for (double ntype : ntypes)
    {
        coverage = ntype/static_cast<double>(nsite);
        coverages.push_back(coverage);
    }

    return coverages;
}

