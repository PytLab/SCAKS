/*******************************************************************************
  * Filename  : test_kmc_functions.c

  * Description: testcase for kmc functions.

  * Author: shaozhengjiang<shaozhengjiang@gmail.com>

  * CreateDate: 2015.12.30

  * ModifyDate: 2015.12.30
*********************************************************************************/

#ifndef TEST_KMC_FUNCTIONS_
#include "test_kmc_functions.h"
#endif

#define NTESTS 3

/* test functions */
void test_collect_coverage(void)
{
    char * possible_types[3] = {"CO", "O", "Vac"};
    double * cvgs;

    cvgs = (double *)malloc(3*sizeof(double));
    collect_coverage(types, 100, possible_types, 3, cvgs, 3);

    // asserts
    CU_ASSERT_DOUBLE_EQUAL(cvgs[0], 0.71, TOLERANCE);
    CU_ASSERT_DOUBLE_EQUAL(cvgs[1], 0.26, TOLERANCE);
    CU_ASSERT_DOUBLE_EQUAL(cvgs[2], 0.03, TOLERANCE);

    free(cvgs);
}


void test_match_elements_list(void)
{
    char * elements_list[4] = {"CO", "O", "CO", "O"};
    double coordinates_list[2][2][3] = {
        {{0.0, 0.0, 0.0},
         {0.0, 1.0, 0.0}},
        {{0.0, 0.0, 0.0},
         {0.0, 1.0, 0.0}}
    };
    int n_success, i;

    n_success = match_elements_list(types, 2, 2, elements_list,
                                    2, 2, 3,
                                    coordinates_list, grid_shape);

    CU_ASSERT_EQUAL(n_success, 30);
}


void test_match_elements(void)
{
    char * elements[2] = {"CO", "O"};
    double coordinates[2][3] = {
        {0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}
    };
    int n_success;

    n_success = match_elements(types, elements, 2, 3, coordinates, grid_shape);

    CU_ASSERT_EQUAL(n_success, 15);
}


/* 
  The suite initialization function.
  returns 0.
*/
int init_suite(void)
{
    return 0;
}


/*
  The suite cleanup function.
  returns 0.
*/
int clean_suite(void)
{
    return 0;
}


/* test case array */
My_TestInfo tests[NTESTS] = {
    {"test collect_coverage", test_collect_coverage},
    {"test match elements", test_match_elements},
    {"test match_elements_list", test_match_elements_list}
};


int main(void)
{
    int i;                    // loop counter
    CU_pSuite pSuite = NULL;  // suit struct pointer

    // initialize the CUnit test registry
    if(CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    // add a suite to the registry
    pSuite = CU_add_suite("kmc_functions_suite", init_suite, clean_suite);
    if(NULL == pSuite)
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    // add the test to suite
    for(i = 0; i < NTESTS; ++i)
    {
        if((NULL == CU_add_test(pSuite, (tests+i)->test_name, (tests+i)->pfunc)))
        {
            CU_cleanup_registry();
            return CU_get_error();
        }
    }

    // run all tests using CUnit
    /*Automated mode*/
    CU_set_output_filename("Test_kmc_functions");
    CU_list_tests_to_file();
    CU_automated_run_tests();

    /*Basic mode*/
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();

    /*Console mode
    CU_console_run_tests();
    */

    CU_cleanup_registry();

    return CU_get_error();
}
