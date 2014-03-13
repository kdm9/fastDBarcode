/*
 * ============================================================================
 *
 *       Filename:  test.c
 *
 *    Description:  tests of fdb
 *
 *        Version:  1.0
 *        Created:  07/03/14 23:05:35
 *       Revision:  none
 *        License:  GPLv3+
 *       Compiler:  gcc
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au [include word penguin in subject]
 *
 * ============================================================================
 */
#include <stdlib.h>
#include "tinytest.h"
#include "tinytest_macros.h"


struct testcase_t fdb_tests[] = {
    END_OF_TESTCASES
};

struct testgroup_t test_groups[] = {
    {"fdb_internals", fdb_tests},
    END_OF_GROUPS
};


int
main (int argc, const char *argv[])
{
    return tinytest_main(argc, argv, test_groups);
} /* ----------  end of function main  ---------- */
