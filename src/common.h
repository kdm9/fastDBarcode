/*
 * ============================================================================
 *
 *       Filename:  common.h
 *
 *    Description:  Some common defs or functions.
 *
 *        Version:  1.0
 *        Created:  15/12/13 22:18:58
 *       Revision:  none
 *        License:  GPLv3+
 *       Compiler:  gcc
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */
#ifndef COMMON_H
#define COMMON_H

#include <errno.h>
#include <libgen.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#define BREAK_EVERY_X_SEQS 1000000

#define FDB_FP_TYPE gzFile
#define FDB_FP_READ gzread
#define FDB_FP_OPEN gzopen
#define FDB_FP_CLOSE gzclose
#define FDB_FP_WRITE gzwrite
#define FDB_FP_ZIP_EXT "gz"

#define FDB_MEM_ERROR \
    fprintf(stderr, "Memory Error: Could not allocate memorary at line %i in %s\n", \
            __LINE__, __FILE__);

#define FDB_IO_ERROR(fle) \
    fprintf(stderr, "IO Error: Could not open file '%s' at line %i in %s\n%s\n", \
            fle, __LINE__, __FILE__, strerror(errno)); \
    exit(EXIT_FAILURE);

#define FDB_MEM_CHECK(var) \
        if (var == NULL) { FDB_MEM_ERROR; }

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

#endif /* COMMON_H */
