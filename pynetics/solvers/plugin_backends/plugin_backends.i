%module plugin_backends

%{
    #define  SWGI_FILE_WITH_INIT
    #include "plugin_backends.h"
%}

%include "numpy.i"

%{
    import_array();
%}
