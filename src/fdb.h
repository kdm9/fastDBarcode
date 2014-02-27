/*
 * ============================================================================
 *
 *       Filename:  fdb.h
 *
 *    Description:  FasttDBarcode functions
 *
 *        Version:  1.0
 *        Created:  16/02/14 11:33:26
 *       Revision:  none
 *        License:  GPLv3+
 *       Compiler:  gcc
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */
#ifndef FDB_H
#define FDB_H

#include "common.h"

extern size_t hamming_max(const char *seq1, const char *seq2, size_t max);
barcode_t **parse_barcode_file(char *barcode_file, int *num, int flag);
int setup_input(char **infiles, int n_infiles,
        char *out_dir, char *out_mode, char *leftover_suffix, int flag,
        char **infile_basenames, char **infile_exts,
        char **outfile_dirs, FDB_FP_TYPE *leftover_outfps);
extern int cmp_barcode_t_rev(const void *left, const void *right);
#endif /* FDB_H */
