#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#define STRLEN 10  // length of single element


int match_elements(char ** types, char ** elements,
                   int nlocal, int ncomp, double coordinates[nlocal][ncomp],
                   const int grid_shape[2]);

int match_elements_list(char ** types, int nrow, int ncol,
                        char * elements_list[nrow][ncol],
                        double coordinates_list[nrow][ncol][3],
                        const int * grid_shape);
