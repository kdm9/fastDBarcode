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

#include <errno.h>
#include <libgen.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "kdm.h"

#define BREAK_EVERY_X_SEQS 1000000

#define FDB_FP_TYPE gzFile
#define FDB_FP_READ gzread
#define FDB_FP_OPEN gzopen
#define FDB_FP_CLOSE gzclose
#define FDB_FP_WRITE gzwrite
#define FDB_FP_ZIP_EXT "gz"

#define FDB_VERSION "v0.0.1a"

#define	FLG_VERBOSE 1 << 0
#define	FLG_ZIPPED_OUT 1 << 1
#define	FLG_VERY_VERBOSE 1 << 2

#define FDB_ZIP_MODE "w9"
#define FDB_NONZIP_MODE "wT"

/* Don't enforce same-length needle and haystack hamming distance. */
#define FDB_HAMMING_MODE_FROMSTART

#include "kseq.h"
KSEQ_INIT(FDB_FP_TYPE, FDB_FP_READ)

typedef struct __barcode_t {
    kstring_t name;
    kstring_t seq;
    uint64_t count;
    FDB_FP_TYPE *fps;
    char **fns;
} barcode_t;

typedef struct __fdb_config_t {
    int flag;
    char **infns;
    int n_infs;
    char *out_dir;
    char *out_mode;
    char *leftover_suffix;
    char **infn_bases;
    char **infn_exts;
    char **outf_dirs;
    FDB_FP_TYPE *leftover_outfps;
    kseq_t **in_kseqs;
    char *out_suffix;
    char *barcode_file;
    barcode_t **barcodes;
    size_t n_barcodes;
    size_t n_infiles;
    int max_barcode_mismatches;
    int max_buffer_mismatches;
    char *buffer_seq;
    size_t *reads_processed;
} fdb_config_t;

#define FDB_IO_ERROR(fle) \
    fprintf(stderr, "IO Error: Could not open file '%s' at line %i in %s\n%s\n", \
            fle, __LINE__, __FILE__, strerror(errno));

extern size_t hamming_max (const char *seq1, const char *seq2, size_t max);
extern int cmp_barcode_t_rev (const void *left, const void *right);
int parse_args (fdb_config_t *cfg, int argc, char **argv);
int parse_barcode_file (fdb_config_t *cfg);
int setup_files (fdb_config_t *cfg);
int fdb_main (fdb_config_t *cfg);
int fdb_config_destroy (fdb_config_t *cfg);
int print_usage();

#endif /* FDB_H */
