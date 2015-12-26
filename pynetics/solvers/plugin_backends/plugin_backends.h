#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#define STRLEN 10  // length of single element


int match_elements(char ** types, int n_local,
                   char ** elements, 
                   float (*coordinates)[3],
                   const int * grid_shape);


int match_elements_list(char ** types, int nrow, int ncol,
                        char * elements_list[nrow][ncol],
                        float coordinates_list[nrow][ncol][3],
                        const int * grid_shape);
