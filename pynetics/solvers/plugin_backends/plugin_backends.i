/* -*- C -*-  (not really, but good for syntax highlighting) */
/*
  *SWIG interface code for C extesion generation of TOFAnalysis Plugin.
  
  *Author: shaozhengjiang<shaozhengjiang@gmail.com>
  *CreateDate: 2015.12.27
  *ModifyDate: 2015.12.28
*/

%module plugin_backends

%{
    #define  SWIG_FILE_WITH_INIT
    #include "plugin_backends.h"
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
        for(i = 0; i < size; i++)
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


/* treat int[2] as a tuple of int in python */
%typemap(in) int[2](int temp[2]) {   // temp[4] becomes a local variable
    int i;
    if (PyTuple_Check($input))
    {
        if (!PyArg_ParseTuple($input,"ii",temp,temp+1))
        {
            PyErr_SetString(PyExc_TypeError,"tuple must have 2 elements");
            return NULL;
        }
        $1 = &temp[0];
    }
    else
    {
        PyErr_SetString(PyExc_TypeError,"expected a tuple.");
        return NULL;
    }
}

/* python list of string to 1D string array*/
%typemap(in) char * [ANY](char * temp[$1_dim0]){
    int i;
    if(!PySequence_Check($input))
    {
        PyErr_SetString(PyExc_TypeError, "Expecting a sequence witeh $1_dim0 elements");
        return NULL;
    }
    if (PyObject_Length($input) != $1_dim0)
    {
        PyErr_SetString(PyExc_ValueError,"Expecting a sequence with $1_dim0 elements");
        return NULL;
    }

    int size = PyObject_Length($input);

    for(i = 0; i < size; i++)
    {
        PyObject * o = PySequence_GetItem($input, i);
        if(!PyString_Check(o))
        {
            Py_XDECREF(o);
            PyErr_SetString(PyExc_ValueError,"Expecting a sequence of strings");
            return NULL;
        }
        temp[i] = PyString_AsString(o);
        Py_DECREF(o);
    }
    $1 = &temp[0];
}

/* declaration of C functions*/

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
                    stripped_elements_list,
                    stripped_coordinates_list,
                    grid_shape)

Parameters:
-----------
types: The site types at the lattice points as a list, list of str.

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
