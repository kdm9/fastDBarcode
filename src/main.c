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

#include "common.h"
#include "fdb.h"

static int flag = 0;


/*
 * ===  FUNCTION  =============================================================
 *         Name:  print_usage
 *  Description:  Prints the usage instructions for fastDBarcode
 * ============================================================================
 */
int
print_usage                    ()
{                                                                   /* {{{ */
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
}                                                                   /* }}} */



/*
 * ===  FUNCTION  =============================================================
 *         Name:  main
 * ============================================================================
 */
int
main                           (int             argc,
                                char           *argv[])
{                                                                   /* {{{ */
    int max_barcode_mismatches = 1;
    int max_buffer_mismatches = 0;
    int n_infiles = 0;
    int n_barcodes = 0;
    char **infiles = NULL;
    kseq_t **infile_kseqs = NULL;
    char *buffer_seq = NULL;
    char *out_dir = NULL;
    char *out_suffix = NULL;
    char *barcode_file = NULL;
    char *leftover_suffix = NULL;
    barcode_t **barcodes = NULL;
    size_t reads_processed = 0;
    char *out_mode;
    /* Parse all arguments {{{ */
    char c;
    while ((c = getopt(argc, argv, "hvzm:M:B:s:o:l:")) != -1) {
        switch (c) {
            case 'm':
                max_barcode_mismatches = atoi(optarg);
                break;
            case 'M':
                max_buffer_mismatches = atoi(optarg);
                break;
            case 'B':
                buffer_seq = strdup(optarg);
                break;
            case 's':
                out_suffix = strdup(optarg);
                break;
            case 'o':
                out_dir = strdup(optarg);
                break;
            case 'l':
                leftover_suffix = strdup(optarg);
                break;
            case 'z':
                flag |= FLG_ZIPPED_OUT;
                break;
            case 'v':
                if (! flag & FLG_VERBOSE) {
                    flag |= FLG_VERBOSE;
                }
                else {
                    flag |= FLG_VERY_VERBOSE;
                }
                break;
            case 'h':
                print_usage();
                exit(EXIT_SUCCESS);
                break;
            case '?':
                fprintf(stderr, "Bad argument -%c\n", c);
                print_usage();
                exit(EXIT_FAILURE);
        }
    }
    if (flag & FLG_VERBOSE) {
        printf("Being verbose.\n");
    }
    int arg_index = optind;
    if ((arg_index + 1) < argc) {
        FDB_FP_TYPE infile_ptr;
        barcode_file = strdup(argv[arg_index++]);

        barcodes = parse_barcode_file(barcode_file, &n_barcodes, flag);
        n_infiles = argc - arg_index;

        infiles = calloc(n_infiles, sizeof(*infiles));
        FDB_MEM_CHECK(infiles);

        infile_kseqs = calloc(n_infiles, sizeof(*infile_kseqs));
        FDB_MEM_CHECK(infile_kseqs);


        for (int infile_index = 0; infile_index < n_infiles; infile_index++) {
            infiles[infile_index] = strdup(argv[arg_index++]);
            FDB_MEM_CHECK(infiles[infile_index]);

            infile_ptr = FDB_FP_OPEN(infiles[infile_index], "r");
            if (infile_ptr == NULL) {
                FDB_IO_ERROR(infiles[infile_index]);
            }

            infile_kseqs[infile_index] = kseq_init(infile_ptr);

            if (flag & FLG_VERBOSE) {
                printf("Using '%s' as an input file\n", infiles[infile_index]);
            }
        }
    } else {
        fprintf(stderr, "ERROR: insufficent number of arguments\n");
        print_usage();
        return EXIT_FAILURE;
    }
    if (leftover_suffix == NULL) {
        leftover_suffix = strdup("_leftover");
    }
    out_mode = (flag & FLG_ZIPPED_OUT)? FDB_ZIP_MODE: FDB_NONZIP_MODE;
    /* End of argument parsing }}} */

    /* Setup input files {{{ */
    char **infile_basenames = calloc(n_infiles, sizeof(*infile_basenames));
    FDB_MEM_CHECK(infile_basenames);
    char **infile_exts = calloc(n_infiles, sizeof(*infile_exts));
    FDB_MEM_CHECK(infile_exts);
    char **outfile_dirs = calloc(n_infiles, sizeof(*outfile_dirs));
    FDB_MEM_CHECK(outfile_dirs);
    FDB_FP_TYPE *leftover_outfps = calloc(n_infiles, sizeof(*leftover_outfps));
    FDB_MEM_CHECK(leftover_outfps);

    if (setup_input(infiles, n_infiles, out_dir, out_mode, leftover_suffix, flag,
            infile_basenames, infile_exts, outfile_dirs, leftover_outfps)) {
        fprintf(stderr, "[main] ERROR: could not setup input\n");
        goto clean;
    }

    /* Setup output files {{{ */
    for (int bbb = 0; bbb < n_barcodes; bbb++) {
        barcodes[bbb]->fps = calloc(n_infiles, sizeof(*barcodes[bbb]->fps));
        barcodes[bbb]->fns = calloc(n_infiles, sizeof(*barcodes[bbb]->fns));
        for (int fff = 0; fff < n_infiles; fff++) {
            char temp2[2<<16];
            if (infile_exts[fff] != NULL) {
                snprintf(temp2, (2<<16), "%s/%s_%s.%s\0", outfile_dirs[fff],
                        infile_basenames[fff], barcodes[bbb]->name.s,
                        infile_exts[fff]);
            } else {
                snprintf(temp2, (2<<16), "%s/%s_%s\0", outfile_dirs[fff],
                        infile_basenames[fff], barcodes[bbb]->name.s);
            }
            barcodes[bbb]->fns[fff] = strdup(temp2);
            barcodes[bbb]->fps[fff] = FDB_FP_OPEN(barcodes[bbb]->fns[fff],
                    out_mode);
            if (barcodes[bbb]->fps[fff] == NULL) {
                fprintf(stderr, "ERROR: Could not open output file '%s'\n",
                        temp2);
                exit(EXIT_FAILURE);
            }
            if (flag & FLG_VERY_VERBOSE) {
                printf("outfile for %s with barcode %s is %s (bcd #%i)\n",
                        barcodes[bbb]->fns[fff], barcodes[bbb]->name.s,
                        barcodes[bbb]->fps[fff], bbb);
            }
        }
    } /* End of setup of output files }}} */

    /* Main Loop: for each file, split by barcode and write {{{ */
    for (int fff = 0; fff < n_infiles; fff++) {
        printf("Processing %s:\t", infiles[fff]); fflush(stdout);
        reads_processed = 0;
        kseq_t *seq = infile_kseqs[fff];
        while(kseq_read(seq) >= 0) {
            reads_processed++;
            size_t *scores = calloc(n_barcodes, sizeof(*scores));
            char *out_seq = NULL;
            int out_len = 0;
            size_t best_score = SIZE_MAX;
            int best_bcd = 0;
            int best_bcd_len = 0;
            int buffer_match = 0;

            for (int bbb = 0; bbb < n_barcodes; bbb++) {
                barcode_t *bcd = barcodes[bbb];
                scores[bbb] = hamming_max(bcd->seq.s, seq->seq.s,
                        max_barcode_mismatches + 1);
            }
            for (int bbb = 0; bbb < n_barcodes; bbb++){
                if (buffer_seq != NULL) {
                    size_t buffer_hamdist = hamming_max(buffer_seq,
                            seq->seq.s + barcodes[bbb]->seq.l,
                            max_buffer_mismatches + 1);
                    buffer_match = buffer_hamdist <= max_buffer_mismatches;
                } else {
                    /* if no buffer seq, always match */
                    buffer_match = 1;
                }
                if (scores[bbb] <= best_score && \
                        barcodes[bbb]->seq.l >= best_bcd_len && \
                        buffer_match) {
                    best_bcd = bbb;
                    best_bcd_len = barcodes[bbb]->seq.l;
                    best_score = scores[bbb];
                }
            }
            FDB_FP_TYPE this_out_fp;
            if (best_score < max_barcode_mismatches) {
                this_out_fp = barcodes[best_bcd]->fps[fff];
                barcodes[best_bcd]->count++;
            } else {
                best_bcd_len = 0;
                this_out_fp = leftover_outfps[fff];
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
            if (flag & FLG_VERY_VERBOSE) {
                if (best_score < max_barcode_mismatches) {
                    printf("seq %s is from barcode %s with score of %zu.\n",
                            seq->name.s, barcodes[best_bcd]->name.s,
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
            if (reads_processed % BREAK_EVERY_X_SEQS == 0) {
                printf("."); fflush(stdout);
                qsort((void*)barcodes, n_barcodes, sizeof(*barcodes), cmp_barcode_t_rev);
            }
        }
        printf(" done!\n");
        if (flag & FLG_VERBOSE) {
            printf("Processed %zu sequences from %s\n", reads_processed,
                    infiles[fff]);
        }
    } /*  End of main loop }}} */

    if (flag & FLG_VERBOSE) {
        printf("\n\n------------------------------------------------\n");
        printf("[main] Summary of barcodes (reads from all input files):\n");
        for (int ccc = 0; ccc<n_barcodes; ccc++) {
            printf("%s: %d\n", barcodes[ccc]->name.s, barcodes[ccc]->count);
        }
    }
    /* Clean up everything that has been alloc'd {{{ */
clean:
    if (infiles != NULL){
        for (int iii = 0; iii < n_infiles; iii++) {
            if (infiles[iii] != NULL) free(infiles[iii]);
        }
        free(infiles);
    }
    if (infile_kseqs != NULL){
        for (int iii = 0; iii < n_infiles; iii++) {
            if (infile_kseqs[iii] != NULL) {
                FDB_FP_CLOSE(infile_kseqs[iii]->f->f);
                kseq_destroy(infile_kseqs[iii]);
            }
        }
        free(infile_kseqs);
    }
    if (barcode_file != NULL) free(barcode_file);
    if (buffer_seq != NULL) free(buffer_seq);
    if (out_dir != NULL) free(out_dir);
    if (out_suffix != NULL) free(out_suffix);
    if (leftover_suffix != NULL) free(leftover_suffix);
    if (barcodes != NULL) {
        for (int iii = 0; iii < n_barcodes; iii++) {
            if (barcodes[iii] != NULL) {
                if (barcodes[iii]->name.s != NULL) {
                    free(barcodes[iii]->name.s);
                }
                if (barcodes[iii]->seq.s != NULL) {
                    free(barcodes[iii]->seq.s);
                }
                for (int jjj = 0; jjj < n_infiles; jjj++) {
                    if (barcodes[iii]->fns[jjj] != NULL) {
                        free(barcodes[iii]->fns[jjj]);
                    }
                    if (barcodes[iii]->fps[jjj] != NULL) {
                        FDB_FP_CLOSE(barcodes[iii]->fps[jjj]);
                    }
                }
            free(barcodes[iii]);
            }
        }
        free(barcodes);
    }
    if (leftover_outfps != NULL) {
        for (int iii = 0; iii <  n_infiles; iii++) {
            if (leftover_outfps[iii] != NULL) {
                FDB_FP_CLOSE(leftover_outfps[iii]);
            }
        }
        free(leftover_outfps);
    }
    /* End of clean-up code }}}*/

    return EXIT_SUCCESS;
} /* ----------  end of function main  ---------- }}} */
