#ifndef PLUGIN_BACKENDS_
    #define PLUGIN_BACKENDS_
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#define STRLEN 10       // length of single element

extern int match_elements(char ** types, char ** elements,
                          int nlocal, int ncomp, double coordinates[nlocal][ncomp],
                          const int grid_shape[2]);

extern int match_elements_list(char ** types, int nrow, int ncol,
                               char ** elements_list,
                               int dim0, int dim1, int dim2,
                               double coordinates_list[dim0][dim1][dim2],
                               const int grid_shape[2]);

extern void collect_coverage(char ** types, int ntot, char ** possible_types,
                             int nsp, double * cvgs, int ncvgs);
