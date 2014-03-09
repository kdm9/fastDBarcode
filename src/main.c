/*
 * ============================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  A fast barcode splitter
 *
 *        Version:  1.0
 *        Created:  14/12/13 23:18:22
 *       Revision:  none
 *        License:  GPLv3+
 *       Compiler:  gcc
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include "fdb.h"

/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 * ============================================================================
 */
int
main (int argc, char *argv[])
{
    fdb_config_t *cfg = km_calloc(1, sizeof(*cfg), &km_onerr_print);
    /* Parse all arguments */
    if (!parse_args(cfg, argc, argv)) {
        fprintf(stderr, "[main] ERROR: could not parse arguments\n");
        fdb_config_destroy(cfg);
        return EXIT_FAILURE;
    }
    /* Parse barcode file */
    if (!parse_barcode_file(cfg)) {
        fprintf(stderr, "[main] ERROR: could not parse barcode file\n");
        fdb_config_destroy(cfg);
        return EXIT_FAILURE;
    }
    /* Setup input/out files */
    if (!setup_files(cfg)) {
        fprintf(stderr, "[main] ERROR: could not setup input\n");
        fdb_config_destroy(cfg);
        return EXIT_FAILURE;
    }
    /* Setup input/out files */
    if (!fdb_main(cfg)) {
        fprintf(stderr, "[main] ERROR: processing files failed\n");
        fdb_config_destroy(cfg);
        return EXIT_FAILURE;
    }
    fdb_config_destroy(cfg);
    return EXIT_SUCCESS;
} /* ----------  end of function main  ---------- */
