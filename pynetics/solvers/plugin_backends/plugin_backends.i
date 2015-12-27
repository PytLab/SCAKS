/* -*- C -*-  (not really, but good for syntax highlighting) */
/*
  *SWIG interface code for C extesion generation of TOFAnalysis Plugin.
  
  *Author: shaozhengjiang<shaozhengjiang@gmail.com>
  *Date  : 2015.12.27
*/

%module plugin_backends

%{
    #define  SWGI_FILE_WITH_INIT
    #include "plugin_backends.h"
%}

%include "numpy.i"

%{
    import_array();
%}

/* customized typemaps*/

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
