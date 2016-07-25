#ifndef PLUGIN_BACKENDS_
#define PLUGIN_BACKENDS_
#endif

#include <vector>
#include <string>

// STL container iterator alias.
typedef std::vector<double>::iterator DoubleIterType;
typedef std::vector<std::string>::iterator StrIterType;
typedef std::vector<std::string>::const_iterator ConstStrIterType;

/*****************************************************************************
  * Function   : collect_coverage

  * Description: collect statistic about species coverages on a grid.

  * Called By  : CoveragesAnalysis.setup(),
                 CoveragesAnalysis.registerStep()

  * Input:
        @types: The site types at the lattice points
        @possible_types: All possible species type
        @coverage_ratios: The coverages ratio for all basis sites.

  * Return:
        @coverages: The coverages for each species in possible_types.

  * Author: shaozhengjiang<shaozhengjiang@gmail.com>
  * Date  : 2015.12.29
  * Update: 2016.07.25
******************************************************************************/

std::vector<double> collect_coverages(const std::vector<std::string> & types,
                                      const std::vector<std::string> & possible_types,
                                      const std::vector<double> & coverage_ratios);

