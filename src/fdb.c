/*
 * ============================================================================
 *
 *       Filename:  fdb.c
 *
 *    Description:  FasttDBarcode functions
 *
 *        Version:  1.0
 *        Created:  16/02/14 11:32:37
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
 *         Name:    hamming_max
 *  Description:    Calculates the hamming distance up until max. Useful for
 *                      non-exact substring finding.
 *                  If mode is not HAMMING_FROMSTART, arguments must be of
 *                      same length. Preprocessor deals with this.
 * Return Value:    size_t: max if max >= hamming dist, else hamming dist
 * ============================================================================
 */
inline size_t
hamming_max                    (const char     *seq1,
                                const char     *seq2,
                                size_t          max)
{                                                                   /* {{{ */
    size_t len = strlen(seq1);
    /* check seq lengths */
#ifndef FDB_HAMMING_MODE_FROMSTART
/* Don't enforce same-length needle and haystack hamming distance. */
    if (len != strlen(seq2)) {
        return(SIZE_MAX);
    }
#endif
    if (strncmp(seq2, seq1, len) == 0){
        /* exact match */
        return 0;
    } else {
        size_t mismatches = 0;
        size_t iii = 0;
        while(iii < len && mismatches < max) {
            if (seq2[iii] != seq1[iii]) {
                mismatches++;
            }
            iii++;
        }
        return mismatches;
    }
}                                                                   /* }}} */

/*
 * ===  FUNCTION  =============================================================
 *         Name:  parse_barcode_file
 *  Description:  Parses a fasta file containing barcode sequences
 * Return Value:  barcode_t **:
 *                An array of pointers to barcode objects
 * ============================================================================
 */
barcode_t **
parse_barcode_file             (char           *barcode_file,
                                int            *num,
                                int             flag)
{                                                                   /* {{{ */
    size_t n_barcodes = 0;
    size_t alloced_barcodes = 2;
    FDB_FP_TYPE fp = NULL;
    kseq_t * ksq = NULL;
    barcode_t **barcodes = calloc(alloced_barcodes, sizeof(*barcodes));
    fp = FDB_FP_OPEN(barcode_file, "r");
    ksq = kseq_init(fp);
    while (kseq_read(ksq) >= 0) {
        if (ksq->seq.l)
        {
            size_t iii = n_barcodes++;
            /* != 0 below is to protect against int overflow */
            if (n_barcodes >= alloced_barcodes && n_barcodes != 0) {
#ifdef FDB_DEBUG
                printf("reallocing barcodes: from %zu to %zu\n",
                        alloced_barcodes, alloced_barcodes <<1);
#endif
                alloced_barcodes = alloced_barcodes << 1;
                barcodes = realloc(barcodes, alloced_barcodes * sizeof(*barcodes));
            }
            if (flag & FLG_VERBOSE) {
                printf("barcode %s is %s\n", ksq->name.s, ksq->seq.s);
            }
            barcode_t *barcode = calloc(1, sizeof(barcode_t));
            barcode->name.l = ksq->name.l;
            barcode->name.m = ksq->name.m;
            barcode->name.s = strdup(ksq->name.s);
            barcode->seq.l = ksq->seq.l;
            barcode->seq.m = ksq->seq.m;
            barcode->seq.s = strdup(ksq->seq.s);
            barcodes[iii] = barcode;
        }
    }
    /* we wan't the multiplication to wrap around to 0 below, as if it does
     * length is > 2^64-1. */
    barcodes = realloc(barcodes, n_barcodes*sizeof(*barcodes));
    *num = n_barcodes;
    kseq_destroy(ksq);
    FDB_FP_CLOSE(fp);
    if (flag & FLG_VERBOSE) {
        printf("Parsed %zu barcodes from %s\n", n_barcodes, barcode_file);
    }
    return barcodes;
} /* -----  end of function parse_barcode_file  ----- }}} */


int
setup_input                    (char          **infiles,
                                int             n_infiles,
                                char           *out_dir,
                                char           *out_mode,
                                char           *leftover_suffix,
                                int             flag,
                                char          **infile_basenames,
                                char          **infile_exts,
                                char          **outfile_dirs,
                                FDB_FP_TYPE    *leftover_outfps)
{
    for (int fff = 0; fff < n_infiles; fff++) {
        /* base/dirname have to work on a copy of str, it gets mangled*/
        char *infile = strdup(infiles[fff]);
        char *infile_base = strdup(basename(infile));
        char *infile_dir = strdup(dirname(infile));
        char *infile_ext = NULL;
        char *temp = NULL;
        /* restore infile to be a copy of the current input file */
        free(infile);
        infile = strdup(infiles[fff]);
        infile_base = basename(infile_base);
        temp = strchr(infile_base, '.');
        if (temp != NULL) {
            int ext_offset = temp - infile_base;
            temp = strdup(infile_base);
            temp[ext_offset] = '\0';
            free(infile_base);
            infile_base = strdup(temp);
            infile_ext = strdup(temp + ext_offset + 1);
            free(temp);
        }
        if (infile_ext != NULL && flag & FLG_ZIPPED_OUT) {
            infile_ext = strcat(infile_ext, ".");
            infile_ext = strcat(infile_ext, FDB_FP_ZIP_EXT);
        }
        if (infile_ext == NULL) {
            infile_ext = "";
        }
#ifdef  FDB_DEBUG
        printf("the basename of %s is %s\next = %s\n",
               infile, infile_base, infile_ext);
#endif
        if (out_dir == NULL) {
            out_dir = strdup(infile_dir);
        }
        infile_basenames[fff] = infile_base;
        infile_exts[fff] = infile_ext;
        outfile_dirs[fff] = out_dir;
        /* 3 = number of slashes/dots, + 1 \0 */
        size_t leftover_name_len = strlen(out_dir) + strlen(infile_base) + \
                   strlen(leftover_suffix) + strlen(infile_ext) + 3 + 1;
        temp = calloc(leftover_name_len, sizeof(*temp));
        snprintf(temp, leftover_name_len - 1, "%s/%s%s.%s", out_dir,
                infile_base, leftover_suffix, infile_ext);
        leftover_outfps[fff] = FDB_FP_OPEN(temp, out_mode);
        free(temp);
        free(infile);
    }
    return 0;
}

inline int
cmp_barcode_t_rev(const void *left, const void *right)
{
    barcode_t *bcd_l = *((barcode_t **)left), *bcd_r = *((barcode_t**)right);
    return ((int)bcd_r->count - (int)bcd_l->count);
}
