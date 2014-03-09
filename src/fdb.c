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
/* Don't enforce same-length needle and haystack hamming distance. */
#ifndef FDB_HAMMING_MODE_FROMSTART
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
}


/*
 * ===  FUNCTION  =============================================================
 *         Name:  parse_barcode_file
 *  Description:  Parses a fasta file containing barcode sequences
 * Return Value:  int: 1 on success, 0 on failure
 * ============================================================================
 */
int
parse_barcode_file (fdb_config_t *cfg)
{
    size_t alloced_barcodes = 2;
    FDB_FP_TYPE fp = NULL;
    kseq_t * ksq = NULL;
    cfg->barcodes = calloc(alloced_barcodes, sizeof(*(cfg->barcodes)));
    fp = FDB_FP_OPEN(cfg->barcode_file, "r");
    ksq = kseq_init(fp);
    while (kseq_read(ksq) >= 0) {
        if (ksq->seq.l)
        {
            size_t iii = cfg->n_barcodes++;
            /* != 0 below is to protect against int overflow */
            if (cfg->n_barcodes >= alloced_barcodes && cfg->n_barcodes != 0) {
                alloced_barcodes = alloced_barcodes << 1;
                cfg->barcodes = realloc(cfg->barcodes,
                        alloced_barcodes * sizeof(*(cfg->barcodes)));
            }
            if (cfg->flag & FLG_VERBOSE) {
                printf("barcode %s is %s\n", ksq->name.s, ksq->seq.s);
            }
            barcode_t *barcode = calloc(1, sizeof(barcode_t));
            barcode->name.l = ksq->name.l;
            barcode->name.m = ksq->name.m;
            barcode->name.s = strdup(ksq->name.s);
            barcode->seq.l = ksq->seq.l;
            barcode->seq.m = ksq->seq.m;
            barcode->seq.s = strdup(ksq->seq.s);
            cfg->barcodes[iii] = barcode;
        }
    }
    /* we wan't the multiplication to wrap around to 0 below, as if it does
     * length is > 2^64-1. */
    cfg->barcodes = realloc(cfg->barcodes,
            cfg->n_barcodes * sizeof(*(cfg->barcodes)));
    kseq_destroy(ksq);
    FDB_FP_CLOSE(fp);
    if (cfg->flag & FLG_VERBOSE) {
        printf("Parsed %zu barcodes from %s\n",
                cfg->n_barcodes, cfg->barcode_file);
    }
    return 0;
} /* -----  end of function parse_barcode_file  ----- */

/*
 * ===  FUNCTION  =============================================================
 *         Name:  print_usage
 *  Description:  Prints the usage instructions for fastDBarcode
 * ============================================================================
 */
int
print_usage ()
{
    printf("fastDBarcode %s\n\n", FDB_VERSION);
    printf("USAGE:\n");
    printf("\tfastDBarcode [-m -M -B -v -o -s -z] <barcode_file> <fq_file> ...\n\n");
    printf("\tfastDBarcode -h\n\n");
    printf("OPTIONS:\n");
    printf("\t-m BCD_MISMATCH\tThe maximal hamming distance between barcode\n");
    printf("\t\t\tand sequences. [DEFAULT 1]\n");
    printf("\t-M BFR_MISMATCH\tThe hamming distance between post-barcode\n");
    printf("\t\t\tbuffer seq and sequences. [DEFAULT 0]\n");
    printf("\t-B BUFFER_SEQ\tSequence after the barcode to match.\n");
    printf("\t-s\t\tOutfile suffix. [DEFAULT barcode_id]\n");
    printf("\t-l\t\tLeftover file suffix. [DEFAULT \"_leftover\"]\n");
    printf("\t-o\t\tOutput directory. [DEFAULT dirname(input) for each file]\n");
    printf("\t-z\t\tWrite output fastqs as zipped files.\n");
    printf("\t-v\t\tBe more verbose.\n");
    printf("\t-h\t\tProvide some help.\n");
    return EXIT_SUCCESS;
}

/*
 * ===  FUNCTION  =============================================================
 *         Name:  setup_files
 *  Description:  Opens and sets up all input and output files, names, pointers
 * ============================================================================
 */
int
setup_files (fdb_config_t *cfg)
{
    for (int fff = 0; fff < cfg->n_infs; fff++) {
        /* base/dirname have to work on a copy of str, it gets mangled*/
        char *infile = strdup(cfg->infns[fff]);
        char *infile_base = strdup(basename(infile));
        char *infile_dir = strdup(dirname(infile));
        char *infile_ext = NULL;
        char *temp = NULL;
        char *out_dir = NULL;
        /* restore infile to be a copy of the current input file */
        free(infile);
        infile = strdup(cfg->infns[fff]);
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
        if (infile_ext != NULL && cfg->flag & FLG_ZIPPED_OUT) {
            infile_ext = strcat(infile_ext, ".");
            infile_ext = strcat(infile_ext, FDB_FP_ZIP_EXT);
        }
        if (infile_ext == NULL) {
            infile_ext = "";
        }
        if (cfg->out_dir == NULL) {
            out_dir = strdup(infile_dir);
        } else {
            out_dir = cfg->out_dir;
        }
        cfg->infn_bases[fff] = infile_base;
        cfg->infn_exts[fff] = infile_ext;
        cfg->outf_dirs[fff] = out_dir;
        /* 3 = number of slashes/dots, + 1 \0 */
        size_t leftover_name_len = strlen(cfg->out_dir) + strlen(infile_base) + \
                   strlen(cfg->leftover_suffix) + strlen(infile_ext) + 3 + 1;
        temp = calloc(leftover_name_len, sizeof(*temp));
        snprintf(temp, leftover_name_len - 1, "%s/%s%s.%s", out_dir,
                infile_base, cfg->leftover_suffix, infile_ext);
        cfg->leftover_outfps[fff] = FDB_FP_OPEN(temp, cfg->out_mode);
        free(temp);
        free(infile);
    }
    /* Setup output files */
    for (int bbb = 0; bbb < cfg->n_barcodes; bbb++) {
        cfg->barcodes[bbb]->fps = calloc(cfg->n_infs,
                sizeof(*(cfg->barcodes[bbb]->fps)));
        cfg->barcodes[bbb]->fns = calloc(cfg->n_infs,
                sizeof(*(cfg->barcodes[bbb]->fns)));
        for (int fff = 0; fff < cfg->n_infs; fff++) {
            char temp2[2<<16];
            if (cfg->infn_exts[fff] != NULL) {
                snprintf(temp2, (2<<16), "%s/%s_%s.%s\0", cfg->outf_dirs[fff],
                        cfg->infn_bases[fff], cfg->barcodes[bbb]->name.s,
                        cfg->infn_exts[fff]);
            } else {
                snprintf(temp2, (2<<16), "%s/%s_%s\0", cfg->outf_dirs[fff],
                        cfg->infn_bases[fff], cfg->barcodes[bbb]->name.s);
            }
            cfg->barcodes[bbb]->fns[fff] = strdup(temp2);
            cfg->barcodes[bbb]->fps[fff] = FDB_FP_OPEN(cfg->barcodes[bbb]->fns[fff],
                    cfg->out_mode);
            if (cfg->barcodes[bbb]->fps[fff] == NULL) {
                fprintf(stderr, "ERROR: Could not open output file '%s'\n",
                        temp2);
                return EXIT_FAILURE;
            }
            if (cfg->flag & FLG_VERY_VERBOSE) {
                printf("outfile for %s with barcode %s is %s (bcd #%i)\n",
                        cfg->barcodes[bbb]->fns[fff], cfg->barcodes[bbb]->name.s,
                        cfg->barcodes[bbb]->fps[fff], bbb);
            }
        }
    } /* End of setup of output files }}} */
    return 0;
}

inline int
cmp_barcode_t_rev (const void *left, const void *right)
{
    barcode_t *bcd_l = *((barcode_t **)left), *bcd_r = *((barcode_t**)right);
    return ((int)bcd_r->count - (int)bcd_l->count);
}

int
parse_args (fdb_config_t *cfg, int argc, char **argv)
{
    char c;
    while ((c = getopt(argc, argv, "hvzm:M:B:s:o:l:")) != -1) {
        switch (c) {
            case 'm':
                cfg->max_barcode_mismatches = atoi(optarg);
                break;
            case 'M':
                cfg->max_buffer_mismatches = atoi(optarg);
                break;
            case 'B':
                cfg->buffer_seq = strdup(optarg);
                break;
            case 's':
                cfg->out_suffix = strdup(optarg);
                break;
            case 'o':
                cfg->out_dir = strdup(optarg);
                break;
            case 'l':
                cfg->leftover_suffix = strdup(optarg);
                break;
            case 'z':
                cfg->flag |= FLG_ZIPPED_OUT;
                break;
            case 'v':
                if (! cfg->flag & FLG_VERBOSE) {
                    cfg->flag |= FLG_VERBOSE;
                }
                else {
                    cfg->flag |= FLG_VERY_VERBOSE;
                }
                break;
            case 'h':
                print_usage();
                return EXIT_SUCCESS;
                break;
            case '?':
                fprintf(stderr, "Bad argument -%c\n", c);
                print_usage();
                return EXIT_FAILURE;
        }
    }
    if (cfg->flag & FLG_VERBOSE) {
        printf("Being verbose.\n");
    }
    int arg_index = optind;
    if ((arg_index + 1) < argc) {
        FDB_FP_TYPE infile_ptr;
        cfg->barcode_file = strdup(argv[arg_index++]);
        cfg->n_infs = argc - arg_index;
        cfg->infns = km_calloc(cfg->n_infs, sizeof(*(cfg->infns)),
                &km_onerr_print);
        cfg->in_kseqs = km_calloc(cfg->n_infs, sizeof(*(cfg->in_kseqs)),
                &km_onerr_print);
        for (int infile_index = 0; infile_index < cfg->n_infs; infile_index++) {
            cfg->infns[infile_index] = strdup(argv[arg_index++]);
            infile_ptr = FDB_FP_OPEN(cfg->infns[infile_index], "r");
            if (infile_ptr == NULL) {
                FDB_IO_ERROR(cfg->infns[infile_index]);
            }
            cfg->in_kseqs[infile_index] = kseq_init(infile_ptr);
            if (cfg->flag & FLG_VERBOSE) {
                printf("Using '%s' as an input file\n", cfg->infns[infile_index]);
            }
        }
    } else {
        fprintf(stderr, "ERROR: insufficent number of arguments\n");
        print_usage();
        return EXIT_FAILURE;
    }
    if (cfg->leftover_suffix == NULL) {
        cfg->leftover_suffix = strdup("_leftover");
    }
    cfg->out_mode = (cfg->flag & FLG_ZIPPED_OUT)? FDB_ZIP_MODE: FDB_NONZIP_MODE;
    /* End of argument parsing }}} */
    cfg->infn_bases = km_calloc(cfg->n_infs, sizeof(*(cfg->infn_bases)),
            &km_onerr_print);
    cfg->infn_exts = km_calloc(cfg->n_infs, sizeof(*(cfg->infn_exts)),
            &km_onerr_print);
    cfg->outf_dirs = km_calloc(cfg->n_infs, sizeof(*(cfg->outf_dirs)),
            &km_onerr_print);
    cfg->leftover_outfps = km_calloc(cfg->n_infs,
            sizeof(*(cfg->leftover_outfps)), &km_onerr_print);
    cfg->reads_processed = km_calloc(cfg->n_infs,
            sizeof(*(cfg->reads_processed)), &km_onerr_print);
}


int
fdb_main (fdb_config_t *cfg)
{
    /* Main Loop: for each file, split by barcode and write {{{ */
    for (int fff = 0; fff < cfg->n_infs; fff++) {
        printf("Processing %s:\t", cfg->infns[fff]); fflush(stdout);
        kseq_t *seq = cfg->in_kseqs[fff];
        while(kseq_read(seq) >= 0) {
            cfg->reads_processed[fff]++;
            size_t *scores = calloc(cfg->n_barcodes, sizeof(*scores));
            char *out_seq = NULL;
            int out_len = 0;
            size_t best_score = SIZE_MAX;
            int best_bcd = 0;
            int best_bcd_len = 0;
            int buffer_match = 0;
            for (int bbb = 0; bbb < cfg->n_barcodes; bbb++) {
                barcode_t *bcd = cfg->barcodes[bbb];
                scores[bbb] = hamming_max(bcd->seq.s, seq->seq.s,
                        cfg->max_barcode_mismatches + 1);
            }
            for (int bbb = 0; bbb < cfg->n_barcodes; bbb++){
                if (cfg->buffer_seq != NULL) {
                    size_t buffer_hamdist = hamming_max(cfg->buffer_seq,
                            seq->seq.s + cfg->barcodes[bbb]->seq.l,
                            cfg->max_buffer_mismatches + 1);
                    buffer_match = buffer_hamdist <= cfg->max_buffer_mismatches;
                } else {
                    /* if no buffer seq, always match */
                    buffer_match = 1;
                }
                if (scores[bbb] <= best_score && \
                        cfg->barcodes[bbb]->seq.l >= best_bcd_len && \
                        buffer_match) {
                    best_bcd = bbb;
                    best_bcd_len = cfg->barcodes[bbb]->seq.l;
                    best_score = scores[bbb];
                }
            }
            FDB_FP_TYPE this_out_fp;
            if (best_score < cfg->max_barcode_mismatches) {
                this_out_fp = cfg->barcodes[best_bcd]->fps[fff];
                cfg->barcodes[best_bcd]->count++;
            } else {
                best_bcd_len = 0;
                this_out_fp = cfg->leftover_outfps[fff];
            }
            /* Write out the barcode */
            /* Extra length: @ + space + '+' 4*\n + \0 = 8
             * -2 * best_bcd_len because we're removing barcode
             */
            out_len = seq->seq.l + seq->name.l + seq->qual.l + \
                      seq->comment.l + 8 - (2 * best_bcd_len);
            out_seq = calloc(out_len, sizeof(*out_seq));
            //snprintf(out_seq, out_len, "@%s %s\n%s\n+\n%s\n",
            snprintf(out_seq, out_len, "@%s %s\n%s\n+\n%s\n",
                    seq->name.s, seq->comment.s, seq->seq.s + best_bcd_len,
                    seq->qual.s + best_bcd_len);
            /* Be verbose about things if we're aksed to */
            if (cfg->flag & FLG_VERY_VERBOSE) {
                if (best_score < cfg->max_barcode_mismatches) {
                    printf("seq %s is from barcode %s with score of %zu.\n",
                            seq->name.s, cfg->barcodes[best_bcd]->name.s,
                            best_score);
#ifdef FDB_DEBUG
                    printf("%s\n\n", out_seq);
#endif
                } else {
                    printf("seq %s is from none of the barcodes.\n",
                            seq->name.s);
#ifdef FDB_DEBUG
                    printf("%s\n\n", out_seq);
#endif
                }
            }
            /* out_len-1 because we don't want to write the \0 */
            FDB_FP_WRITE(this_out_fp, out_seq, out_len - 1);
            /* Clean up */
            free(scores);
            free(out_seq);
            if (cfg->reads_processed[fff] % BREAK_EVERY_X_SEQS == 0) {
                printf("."); fflush(stdout);
            }
        }
        printf(" done!\n");
        if (cfg->flag & FLG_VERBOSE) {
            printf("Processed %zu sequences from %s\n",
                    cfg->reads_processed[fff], cfg->infns[fff]);
        }
    } /*  End of main loop }}} */
    if (cfg->flag & FLG_VERBOSE) {
        printf("\n\n------------------------------------------------\n");
        printf("[main] Summary of barcodes (reads from all input files):\n");
        for (int ccc = 0; ccc<cfg->n_barcodes; ccc++) {
            printf("%s: %d\n", cfg->barcodes[ccc]->name.s, cfg->barcodes[ccc]->count);
        }
    }
}

int
fdb_config_destroy (fdb_config_t *cfg)
{
    if (cfg->infns != NULL){
        for (int iii = 0; iii < cfg->n_infs; iii++) {
            if (cfg->infns[iii] != NULL) free(cfg->infns[iii]);
        }
        free(cfg->infns);
    }
    if (cfg->in_kseqs != NULL){
        for (int iii = 0; iii < cfg->n_infs; iii++) {
            if (cfg->in_kseqs[iii] != NULL) {
                FDB_FP_CLOSE(cfg->in_kseqs[iii]->f->f);
                kseq_destroy(cfg->in_kseqs[iii]);
            }
        }
        free(cfg->in_kseqs);
    }
    if (cfg->barcode_file != NULL) free(cfg->barcode_file);
    if (cfg->buffer_seq != NULL) free(cfg->buffer_seq);
    if (cfg->out_dir != NULL) free(cfg->out_dir);
    if (cfg->out_suffix != NULL) free(cfg->out_suffix);
    if (cfg->leftover_suffix != NULL) free(cfg->leftover_suffix);
    if (cfg->barcodes != NULL) {
        for (int iii = 0; iii < cfg->n_barcodes; iii++) {
            if (cfg->barcodes[iii] != NULL) {
                if (cfg->barcodes[iii]->name.s != NULL) {
                    free(cfg->barcodes[iii]->name.s);
                }
                if (cfg->barcodes[iii]->seq.s != NULL) {
                    free(cfg->barcodes[iii]->seq.s);
                }
                for (int jjj = 0; jjj < cfg->n_infs; jjj++) {
                    if (cfg->barcodes[iii]->fns[jjj] != NULL) {
                        free(cfg->barcodes[iii]->fns[jjj]);
                    }
                    if (cfg->barcodes[iii]->fps[jjj] != NULL) {
                        FDB_FP_CLOSE(cfg->barcodes[iii]->fps[jjj]);
                    }
                }
            free(cfg->barcodes[iii]);
            }
        }
        free(cfg->barcodes);
    }
    if (cfg->leftover_outfps != NULL) {
        for (int iii = 0; iii <  cfg->n_infs; iii++) {
            if (cfg->leftover_outfps[iii] != NULL) {
                FDB_FP_CLOSE(cfg->leftover_outfps[iii]);
            }
        }
        free(cfg->leftover_outfps);
    }
}
