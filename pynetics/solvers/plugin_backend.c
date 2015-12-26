#include "plugin_backend.h"

int main(void)
{
    char * elements_list[2][2] = {
        {"CO", "O"},
        {"CO", "O"}
    };
    float coordinates_list[2][2][3] = {
        {{0.0, 0.0, 0.0},
         {0.0, 1.0, 0.0}},
        {{0.0, 0.0, 0.0},
         {0.0, 1.0, 0.0}}
    };
    char * types[100] = {
        "O", "O", "CO", "O", "O", "CO", "O", "CO", 
        "CO", "CO", "CO", "CO", "CO", "O", "O", "CO", 
        "O", "CO", "CO", "CO", "CO", "CO", "CO", "CO", 
        "CO", "CO", "CO", "CO", "O", "O", "CO", "CO", 
        "CO", "CO", "CO", "CO", "CO", "CO", "CO", "CO", 
        "CO", "O", "O", "CO", "CO", "O", "O", "CO", 
        "O", "O", "CO", "Vac", "CO", "CO", "CO", "CO", 
        "CO", "CO", "CO", "CO", "CO", "CO", "CO", "CO", 
        "CO", "O", "O", "Vac", "CO", "CO", "O", "CO", 
        "CO", "CO", "CO", "CO", "CO", "CO", "CO", "O", 
        "CO", "CO", "CO", "CO", "O", "O", "CO", "CO", 
        "O", "CO", "CO", "CO", "CO", "O", "O", "CO", 
        "CO", "CO", "O", "Vac"
    };
    int grid_shape[2] = {10, 10};
    int n_success;

    n_success = match_elements_list(types, 2, 2, elements_list, coordinates_list, grid_shape);

    printf("n_success = %d\n", n_success);

    return 0;
}


/**************************************************************
  * Function   : match_elements

  * Description: match local configuration in a grid

  * Called By  : match_elements_list()

  * Input:
        @types      : The site types at the lattice points
        @elements   : elements in local configuration
        @coordinates: relative coordinates of local configuration
        @n_local    : number of sites of local configuration
        @grid_shape : shape of lattice surface grid

  * Return:
        @n_success: number of successful matching

  * Author: shaozhengjiang<shaozhengjiang@gmail.com>

  * Date  : 2015.12.26
***************************************************************/
int match_elements(char ** types, int n_local,
                   char ** elements, 
                   float (*coordinates)[3],
                   const int * grid_shape)
{
    int nrow, ncol;            // grid shape, number rows and columns
    int n_success;             // number of successful matching
    int i, j;                  // counters for row and column
    int ilocal;                // counter for local elements
    bool match_fail;           // matching failure ?
    float x_offset, y_offset;  // position offset vector
    int x, y;                  // offseted coordinates components
    int idx;                   // index of types
    char element[STRLEN];      // element type

    nrow = grid_shape[0];
    ncol = grid_shape[1];

    // go through every site on grid
    for(i = 0, n_success = 0; i < nrow; ++i)
    {
        for(j = 0; j < ncol; ++j)
        {
            match_fail = false;
            for(ilocal = 0; ilocal < n_local; ++ilocal)
            {
                // get offset vector components
                x_offset = coordinates[ilocal][0];
                y_offset = coordinates[ilocal][1];

                // get element type
                strcpy(element, elements[ilocal]);

                // do periodic boundary condition check
                x = i + (int)x_offset;
                y = j + (int)y_offset;
                // check x
                if(x < 0)
                    x = nrow - 1;
                else if(x > nrow - 1)
                    x = 0;
                // check y
                if(y < 0)
                    y = ncol - 1;
                else if(y > ncol - 1)
                    y = 0;
                
                // index in types
                idx = x*ncol + y;

                // compare elements with target element
                if(strcmp(types[idx], element))
                {
                    match_fail = true;
                    break;
                }
            }
            // collect success number
            if(!match_fail)
                n_success++;
        }
    }

    return n_success;
}


/*******************************************************************************
  * Function   : match_elements_list

  * Description: Function to get total matching success number
                 for a list of elements and coordinates

  * Called By  : kmc_plugin.TOFAnalysis.registerStep()

  * Input:
        @types      : The site types at the lattice points
        @grid_shape : shape of lattice surface grid
        @nrow : length of 0th dimesion of elements_list and coordinates_list
        @ncol : length of 1st dimesion of elements_list and coordinates_list
        @elements_list   : a list of elements in local configuration
        @coordinates_list: a list of relative coordinates of local configuration

  * Return:
        @total_success: total number of successful matching

  * Author: shaozhengjiang<shaozhengjiang@gmail.com>

  * Date  : 2015.12.26
*********************************************************************************/
int match_elements_list(char ** types, int nrow, int ncol,
                        char * elements_list[nrow][ncol],
                        float coordinates_list[nrow][ncol][3],
                        const int * grid_shape)
{
    int total_success;  // total number of successful matching
    int n_success;      // number of successful matching for a local config
    int irow;           // row counter
    char ** elements;   // elements in local configuration
    float (*coordinates)[3];

    for(total_success = 0, irow = 0; irow < nrow; ++irow)
    {
        // get elements
        elements = elements_list[irow];

        // get coordinates
        coordinates = coordinates_list[irow];

        // match that local configuration
        n_success = match_elements(types, ncol, elements, coordinates, grid_shape);

        // add to total success number
        total_success += n_success;
    }

    return total_success;
}
