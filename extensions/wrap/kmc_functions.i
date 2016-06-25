/* -*- C -*-  (not really, but good for syntax highlighting) */

/*****************************************************************
 * SWIG interface code for C++ extesion generation of KMC Plugin.
   
 * Author: shaozhengjiang<shaozhengjiang@gmail.com>

 * Date  : 2015.12.27

 * Update: 2016.06.25
 *****************************************************************/

%module kmc_functions

// Preprocessor directives in *_wrap.cxx wrap code.
%{
#include "kmc_functions.h"
%}

// Include SWIG interface files for STL containers.
%include "std_vector.i"
%include "std_string.i"

// Define the templates to be used in Python.
%template(StdVectorString) std::vector<std::string>;
%template(StdVectorDouble) std::vector<double>;

// Functions to be wrapped.
%include "kmc_functions.h"

