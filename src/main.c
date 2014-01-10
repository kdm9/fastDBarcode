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
 *         Author:  Kevin Murray, spam@kdmurray.id.au [include word penguin in subject]
 *
 * ============================================================================
 */

#include "common.h"


static int flag = 0;

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
inline static size_t
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
    printf("\t-o\t\tOutput directory. [DEFAULT dirname(input) for each file]\n");
    printf("\t-z\t\tWrite output fastqs as zipped files.\n");
    printf("\t-v\t\tBe more verbose.\n");
    printf("\t-h\t\tProvide some help.\n");
    return EXIT_SUCCESS;
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
                                int            *num)
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
            barcode->name = calloc(1, sizeof(kstring_t));
            barcode->seq = calloc(1, sizeof(kstring_t));
            barcode->name->l = ksq->name.l;
            barcode->name->m = ksq->name.m;
            barcode->name->s = strdup(ksq->name.s);
            barcode->seq->l = ksq->seq.l;
            barcode->seq->m = ksq->seq.m;
            barcode->seq->s = strdup(ksq->seq.s);
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
    FDB_FP_TYPE *infile_ptrs = NULL;
    kseq_t **infile_kseqs = NULL;
    char *buffer_seq = NULL;
    char *out_dir = NULL;
    char *out_suffix = NULL;
    FDB_FP_TYPE *outfile_ptrs = NULL;
    char **outfiles = NULL;
    char *barcode_file = NULL;
    barcode_t **barcodes = NULL;
    size_t reads_processed = 0;

    /* Parse all arguments {{{ */
    char c;
    while ((c = getopt(argc, argv, "hvzm:M:B:s:o:")) != -1) {
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
            case 'z':
                flag |= FLG_ZIPPED_OUT;
                break;
            case 'v':
                flag |= FLG_VERBOSE;
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
        barcode_file = strdup(argv[arg_index++]);

        barcodes = parse_barcode_file(barcode_file, &n_barcodes);
        n_infiles = argc - arg_index;
        infiles = calloc(n_infiles, sizeof(*infiles));
        FDB_MEM_CHECK(infiles);

        infile_ptrs = calloc(n_infiles, sizeof(*infile_ptrs));
        FDB_MEM_CHECK(infile_ptrs)

        infile_kseqs = calloc(n_infiles, sizeof(*infile_kseqs));
        FDB_MEM_CHECK(infile_kseqs);

        for (int infile_index = 0; infile_index < n_infiles; infile_index++) {
            infiles[infile_index] = strdup(argv[arg_index++]);
            FDB_MEM_CHECK(infiles[infile_index]);

            infile_ptrs[infile_index] = FDB_FP_OPEN(infiles[infile_index], "r");
            if (infile_ptrs[infile_index] == NULL) {
                FDB_IO_ERROR(infiles[infile_index]);
            }

            infile_kseqs[infile_index] = kseq_init(infile_ptrs[infile_index]);

            if (flag & FLG_VERBOSE) {
                printf("Using '%s' as an input file\n", infiles[infile_index]);
            }
        }
    } else {
        fprintf(stderr, "ERROR: insufficent number of arguments\n");
        print_usage();
        return EXIT_FAILURE;
    }
    /* End of argument parsing }}} */

    /* Setup output files {{{ */
    outfiles = calloc(n_infiles * n_barcodes, sizeof(*outfiles));
    outfile_ptrs = calloc(n_infiles * n_barcodes, sizeof(*outfile_ptrs));
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

#ifdef  FDB_DEBUG
        printf("the basename of %s is %s\n", infile, infile_base);
#endif

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

#ifdef  FDB_DEBUG
        printf("the basename of %s is %s\next = %s\n",
               infile, infile_base, infile_ext);
#endif
        if (out_dir == NULL) {
            out_dir = strdup(infile_dir);
        }
        for (int bbb = 0; bbb < n_barcodes; bbb++) {
            int ooo = fff * n_barcodes + bbb;
            char temp2[2<<10] = "";
            char *mode = (flag & FLG_ZIPPED_OUT)? \
                          FDB_ZIP_MODE: FDB_NONZIP_MODE;
            if (infile_ext != NULL) {
                snprintf(temp2, (2<<10)-1, "%s/%s_%s.%s", out_dir, infile_base,
                        barcodes[bbb]->name->s, infile_ext);
            } else {
                snprintf(temp2, (2<<10)-1, "%s/%s_%s", out_dir, infile_base,
                        barcodes[bbb]->name->s);
            }
            outfiles[ooo] = strdup(temp2);
            outfile_ptrs[ooo] =  FDB_FP_OPEN(outfiles[ooo], mode);
            if (flag & FLG_VERBOSE) {
                printf("outfile for %s with barcode %s is %s (bcd #%i)\n",
                        infile, barcodes[bbb]->name->s, outfiles[ooo], ooo);
            }
        }
        free(infile);
        free(infile_base);
        free(infile_dir);
        free(infile_ext);
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
                scores[bbb] = hamming_max(bcd->seq->s, seq->seq.s,
                        max_barcode_mismatches + 1);
            }
            for (int bbb = 0; bbb < n_barcodes; bbb++){
                /* printf("idx %i, bcd %s, score %zu, len %zu\n",
                        bbb, barcodes[bbb]->name->s, scores[bbb],
                        barcodes[bbb]->seq->l); */
                if (buffer_seq != NULL) {
                    size_t buffer_hamdist = hamming_max(buffer_seq,
                            seq->seq.s + barcodes[bbb]->seq->l,
                            max_buffer_mismatches + 1);
                    buffer_match = buffer_hamdist <= max_buffer_mismatches;
                    /* printf("bufseq %s, bdist %zu, match %i\n", buffer_seq,
                            buffer_hamdist, buffer_match); */
                } else {
                    /* if no buffer seq, always match */
                    buffer_match = 1;
                }
                if (scores[bbb] <= best_score && \
                        barcodes[bbb]->seq->l >= best_bcd_len && \
                        buffer_match) {
                    best_bcd = bbb;
                    best_bcd_len = barcodes[bbb]->seq->l;
                    best_score = scores[bbb];
                }
           }
            if (best_score < max_barcode_mismatches) {
               /* Extra length: @ + space + '+' 4*\n + \0 = 8
                 * -2 * best_bcd_len because we're removing barcode
                 */
                out_len = seq->seq.l + seq->name.l + seq->qual.l + \
                          seq->comment.l + 8 - 2 * best_bcd_len;
                out_seq = calloc(out_len, sizeof(*out_seq));
                snprintf(out_seq, out_len, "@%s %s\n%s\n+\n%s\n",
                        seq->name.s, seq->comment.s, seq->seq.s + best_bcd_len,
                        seq->qual.s + best_bcd_len);
                /* out_len-1 because we don't want to write the \0, which it w*/
                FDB_FP_WRITE(outfile_ptrs[best_bcd], out_seq, out_len - 1);

                if (flag & FLG_VERBOSE) {
                    printf("seq %s is from barcode %s with score of %zu\n\n%s\n\n",
                            seq->name.s, barcodes[best_bcd]->name->s, best_score,
                            out_seq);
                }
            }
            free(scores);
            free(out_seq);
            if (reads_processed % 1000000 == 0) { printf("."); fflush(stdout); }
        }
        printf(" done!\n");
        if (flag & FLG_VERBOSE) {
            printf("Processed %zu sequences from %s\n", reads_processed,
                    infiles[fff]);
        }
    } /*  End of main loop }}} */

    /* Clean up everything that has been alloc'd {{{ */
    if (infiles != NULL){
        for (int iii = 0; iii < n_infiles; iii++) {
            if (infiles[iii] != NULL) free(infiles[iii]);
        }
        free(infiles);
    }
    if (infile_ptrs != NULL){
        for (int iii = 0; iii < n_infiles; iii++) {
            if (infile_ptrs[iii] != NULL) FDB_FP_CLOSE(infile_ptrs[iii]);
        }
        free(infile_ptrs);
    }
    if (infile_kseqs != NULL){
        for (int iii = 0; iii < n_infiles; iii++) {
            if (infile_kseqs[iii] != NULL) kseq_destroy(infile_kseqs[iii]);
        }
        free(infile_kseqs);
    }
    if (barcode_file != NULL) free(barcode_file);
    if (buffer_seq != NULL) free(buffer_seq);
    if (out_dir != NULL) free(out_dir);
    if (out_suffix != NULL) free(out_suffix);
    if (barcodes != NULL) {
        for (int iii = 0; iii < n_barcodes; iii++) {
            if (barcodes[iii] != NULL) {
                if (barcodes[iii]->name != NULL) {
                    if (barcodes[iii]->name->s != NULL) {
                        free(barcodes[iii]->name->s);
                    }
                    free(barcodes[iii]->name);
                }
                if (barcodes[iii]->seq != NULL) {
                    if (barcodes[iii]->seq->s != NULL) {
                        free(barcodes[iii]->seq->s);
                    }
                    free(barcodes[iii]->seq);
                }
                free(barcodes[iii]);
            }
        }
        free(barcodes);
    }
    if (outfiles != NULL) {
        for (int iii = 0; iii < n_barcodes * n_infiles; iii++) {
            if (outfiles[iii] != NULL) {
                free(outfiles[iii]);
            }
        }
        free(outfiles);
    }
    if (outfile_ptrs != NULL) {
        for (int iii = 0; iii < n_barcodes * n_infiles; iii++) {
            if (outfile_ptrs[iii] != NULL) {
                FDB_FP_CLOSE(outfile_ptrs[iii]);
            }
        }
        free(outfile_ptrs);
    }
    /* End of clean-up code }}}*/

    return EXIT_SUCCESS;
} /* ----------  end of function main  ---------- }}} */
