#ifndef PLUGIN_BACKENDS_
#include "plugin_backends.h"
#endif
#include <omp.h>

/**************************************************************
  * Function   : collect_coverage

  * Description: collect statistic about species coverages on a grid.

  * Called By  : CoveragesAnalysis.setup(),
                 CoveragesAnalysis.registerStep()

  * Input:
        @types: The site types at the lattice points
        @ntot : total number of site on grid
        @possible_types: all possible species type
        @nsp  : number of possible species type
        @cvgs : result, species coverages on this grid
        @ncvgs: equal to nsp

  * Return:
        None

  * Author: shaozhengjiang<shaozhengjiang@gmail.com>

  * Date  : 2015.12.29
***************************************************************/
void collect_coverage(char ** types, int ntot, char ** possible_types,
                      int nsp, double * cvgs, int ncvgs)
{
    int itype, isp;  // loop counter for types and nspecies
    bool same;       // species matching successful or not

    // initialize coverages
    for(isp = 0; isp < nsp; ++isp)
        cvgs[isp] = 0.0;

    // collect species
    for(itype = 0; itype < ntot; ++itype)
    {
        for(isp = 0; isp < nsp; ++isp)
        {
            same = (strcmp(types[itype], possible_types[isp]) == 0);
            if(same) cvgs[isp] += 1.0;
        }
    }

    // get coverges
    for(isp = 0; isp < nsp; ++isp)
        cvgs[isp] = cvgs[isp]/(float)ntot;
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
int match_elements(char ** types, char ** elements,
                   int nlocal, int ncomp, double coordinates[nlocal][ncomp],
                   const int grid_shape[2])
{
    int nrow, ncol;             // grid shape, number rows and columns
    int n_success;              // number of successful matching
    int i, j;                   // counters for row and column
    int ilocal;                 // counter for local elements
    bool match_fail;            // matching failure ?
    double x_offset, y_offset;  // position offset vector
    int x, y;                   // offseted coordinates components
    int idx;                    // index of types
    char element[STRLEN];       // element type

    nrow = grid_shape[0];
    ncol = grid_shape[1];

    // go through every site on grid
    n_success = 0; 
    for(i = 0; i < nrow; ++i)
    {
        for(j = 0; j < ncol; ++j)
        {
            match_fail = false;
            for(ilocal = 0; ilocal < nlocal; ++ilocal)
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
        @nrow : length of 0th dimesion of elements_list
        @ncol : length of 1st dimesion of elements_list
        @elements_list: a 1D array of elements in local configuration
        @dim0 : length of 0th dimension of coordinates_list
        @dim1 : length of 1st dimension of coordinates_list
        @dim2 : length of 2nd dimension of coordinates_list
        @coordinates_list: a list of relative coordinates of local configuration
        @grid_shape : shape of lattice surface grid

  * Return:
        @total_success: total number of successful matching

  * Author: shaozhengjiang<shaozhengjiang@gmail.com>

  * CreateDate: 2015.12.26

  * ModifyDate: 2015.12.28
*********************************************************************************/
int match_elements_list(char ** types, int nrow, int ncol,
                        char ** elements_list,
                        int dim0, int dim1, int dim2,
                        double coordinates_list[dim0][dim1][dim2],
                        const int grid_shape[2])
{
    int total_success;  // total number of successful matching
    int n_success;      // number of successful matching for a local config
    int irow;           // row counter
    char ** elements;   // elements in local configuration
    int index;          // index number in 1D char * array -- elements_list
    double (*coordinates)[dim2];  // point to a 2D int array

    total_success = 0;
#pragma omp parallel for
    for(irow = 0; irow < nrow; ++irow)
    {
        // get elements
        index = ncol * irow;
        elements = &elements_list[index];

        // get coordinates
        coordinates = coordinates_list[irow];

        // match that local configuration
        n_success = match_elements(types, elements, ncol, dim2,
                                   coordinates, grid_shape);

        // add to total success number
        total_success += n_success;
    }


    return total_success;
}
