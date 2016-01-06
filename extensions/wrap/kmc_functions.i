/* -*- C -*-  (not really, but good for syntax highlighting) */
/*
  *SWIG interface code for C extesion generation of TOFAnalysis Plugin.
  
  *Author: shaozhengjiang<shaozhengjiang@gmail.com>
  *CreateDate: 2015.12.27
  *ModifyDate: 2015.12.28
*/

%define DOCSTR
"C extension module for kMC Analysis key functions."
%enddef

%module(docstring=DOCSTR) kmc_functions

%{
    #define  SWIG_FILE_WITH_INIT
    #ifndef PLUGIN_BACKENDS_
        #include "plugin_backends.h"
    #endif
    #include <omp.h>
%}

%include "numpy.i"

%init %{
    import_array();
%}

/* 
  customized typemaps
*/

/* tell SWIG to treat char ** as a list of strings */
%typemap(in) char ** {
    // check if is a list
    if(PyList_Check($input))
    {
        int size = PyList_Size($input);
        int i = 0;
        $1 = (char **)malloc((size + 1)*sizeof(char *));
        for(i = 0; i < size; ++i)
        {
            PyObject * o = PyList_GetItem($input, i);
            if(PyString_Check(o))
                $1[i] = PyString_AsString(o);
            else
            {
                PyErr_SetString(PyExc_TypeError, "list must contain strings");
                free($1);
                return NULL;
            }
        }
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "not a list");
        return NULL;
    }
}

// clean up the char ** array
%typemap(freearg) char ** {
    free((char *) $1);
}


/* typemaps for collect_coverage function */
/* 
   treat string list as a C string array, 
   just take one input in collect_coverage
*/
%typemap(in, numinputs=1) (char ** TYPES, int NTYPES){
    int i;
    assert(PyList_Check($input));
    $2 = PyList_Size($input);
    $1 = (char **)malloc(($2 + 1)*sizeof(char *));

    for(i = 0; i < $2; ++i)
    {
        PyObject * o = PyList_GetItem($input, i);
        if(PyString_Check(o))
            $1[i] = PyString_AsString(o);
        else
        {
            PyErr_SetString(PyExc_TypeError, "list must contain strings");
            free($1);
            return NULL;
        }
    }
}

// clean up the char ** array
%typemap(freearg) (char ** TYPES, int NTYPES) {
    free((char *) $1);
}


%typemap(in, numinputs=1) (char ** POSSIBLE_TYPES, int NSP){
    int i;
    assert(PyList_Check($input));
    $2 = PyList_Size($input);
    $1 = (char **)malloc(($2 + 1)*sizeof(char *));

    for(i = 0; i < $2; ++i)
    {
        PyObject * o = PyList_GetItem($input, i);
        if(PyString_Check(o))
            $1[i] = PyString_AsString(o);
        else
        {
            PyErr_SetString(PyExc_TypeError, "list must contain strings");
            free($1);
            return NULL;
        }
    }
}

// clean up the char ** array
%typemap(freearg) (char ** POSSIBLE_TYPES, int NSP) {
    free((char *) $1);
}


/* treat int[2] as a tuple of int in python */
%typemap(in) int[2](int temp[2]) {   // temp[4] becomes a local variable
    int i;
    if (PyTuple_Check($input))
    {
        if (!PyArg_ParseTuple($input, "ii", temp, temp+1))
        {
            PyErr_SetString(PyExc_TypeError, "tuple must have 2 elements");
            return NULL;
        }
        $1 = &temp[0];
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "expected a tuple.");
        return NULL;
    }
}

/* declaration of C functions*/

/* TOF Analysis */
%define MATCH_ELEMENTS_DOCSTRING
"Function to go through grid to match elements local configuration.

Python function:
----------------
match_elements(types, stripped_elements, stripped_coordinates, grid_shape)

Parameters:
-----------
types: The site types at the lattice points as a list, list of str.

stripped_elements: stripped elements list(without wildcards),
                   numpy.array of str.

stripped_coordinates: stripped relative coordinates list(without wildcards),
                      2d numpy.array of float.

grid_shape: shape of grid, tuple of int.

Returns:
--------
n_success: number of successful matching, int
"
%enddef

%feature("autodoc", MATCH_ELEMENTS_DOCSTRING);
int match_elements(char ** types, char ** elements,
                   int DIM1, int DIM2, double * IN_ARRAY2,
                   const int grid_shape[2]);

%define MATCH_ELEMENTS_LIST_DOCSTRING
"Function to get total matching success number for,
a list of stripped elements list and coordinates.

Python function prototype:
--------------------------
match_elements_list(types,
                    nrow, ncol,
                    stripped_elements_list,
                    stripped_coordinates_list,
                    grid_shape)

Parameters:
-----------
types: The site types at the lattice points as a list, list of str.

nrow, ncol: numbers of row and column of lattice grid.

stripped_elements_list: a list of stripped_elements, a **1D** string list.

stripped_coordinates_list: a list of stripped coordinates, 2D string list.

grid_shape: shape of grid, tuple of int.

Returns:
--------
total_nsuccess: total number of successful matching, int
"
%enddef

%feature("autodoc", MATCH_ELEMENTS_LIST_DOCSTRING);
int match_elements_list(char ** types, int nrow, int ncol,
                        char ** elements_list,
                        int DIM1, int DIM2, int DIM3,
                        double * IN_ARRAY3,
                        const int grid_shape[2]);


/* Coverage Analysis */
%define COLLECT_COVERAGE_DOCSTRING
"Function to get current coverages of possible types.

Python function prototype:
--------------------------
collect_coverage(types, possible_types, ncvgs)

Parameters:
-----------
types: The site types at the lattice points as a list, list of str.

possible_types: possible species type in grid.

ncvgs: number of possible species type.

Returns:
--------
cvgs: coverages of possible types, numpy.array int
"
%enddef

%feature("autodoc", COLLECT_COVERAGE_DOCSTRING);
void collect_coverage(char ** TYPES, int NTYPES, char ** POSSIBLE_TYPES,
                      int NSP, double * ARGOUT_ARRAY1, int DIM1);
