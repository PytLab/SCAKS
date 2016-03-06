#ifndef TEST_KMC_FUNCTIONS_
#define TEST_KMC_FUNCTIONS_
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <CUnit/Basic.h>
#include <CUnit/Console.h>
#include <CUnit/CUnit.h>
#include <CUnit/TestDB.h>

#ifndef PLUGIN_BACKENDS_
#include "plugin_backends.h"
#endif

#define TOLERANCE 1e-10


struct My_TestInfo{
    char * test_name;
    void (*pfunc)(void);
};
typedef struct My_TestInfo My_TestInfo;

typedef My_TestInfo* My_pTestInfo;


/* Global Varibales: */
char * types[100] = {
    "O" , "O" , "CO", "O"  , "O" , "CO", "O" , "CO",
    "CO", "CO", "CO", "CO" , "CO", "O" , "O" , "CO",
    "O" , "CO", "CO", "CO" , "CO", "CO", "CO", "CO",
    "CO", "CO", "CO", "CO" , "O" , "O" , "CO", "CO",
    "CO", "CO", "CO", "CO" , "CO", "CO", "CO", "CO",
    "CO", "O" , "O" , "CO" , "CO", "O" , "O" , "CO",
    "O" , "O" , "CO", "Vac", "CO", "CO", "CO", "CO",
    "CO", "CO", "CO", "CO" , "CO", "CO", "CO", "CO",
    "CO", "O" , "O" , "Vac", "CO", "CO", "O" , "CO",
    "CO", "CO", "CO", "CO" , "CO", "CO", "CO", "O" ,
    "CO", "CO", "CO", "CO" , "O" , "O" , "CO", "CO",
    "O" , "CO", "CO", "CO" , "CO", "O" , "O" , "CO",
    "CO", "CO", "O" , "Vac"
};

int grid_shape[2] = {10, 10};
