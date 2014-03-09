/*
 * ============================================================================
 *
 *       Filename:  kdm.h
 *
 *    Description:  libkdm: header-only utility to make C more DRY-ish
 *
 *        Version:  1.0
 *        Created:  03/03/14 17:06:23
 *       Revision:  none
 *        License:  GPLv3+
 *       Compiler:  gcc / clang
 *
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#ifndef KDM_H
#define KDM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifndef _WIN32
#include <unistd.h>
#endif


/*
 * Macro helpers from tor
 */

/** Expands to a syntactically valid empty statement.  */
#define STMT_NIL (void)0

/** Expands to a syntactically valid empty statement, explicitly (void)ing its
 * argument. */
#define STMT_VOID(a) while (0) { (void)(a); }

#ifdef __GNUC__
/** STMT_BEGIN and STMT_END are used to wrap blocks inside macros so that
 * the macro can be used as if it were a single C statement. */
#define STMT_BEGIN (void) ({
#define STMT_END })
#elif defined(sun) || defined(__sun__)
#define STMT_BEGIN if (1) {
#define STMT_END } else STMT_NIL
#else
#define STMT_BEGIN do {
#define STMT_END } while (0)
#endif

/*
 * Error handling constants
 */

#define KDM_ERR_ALLOC 1<<0
#define KDM_ERR_FREE 1<<1

static const char *km_err_msgs[] = {
    "No Error",
    "Could not allocate memory",
    "Could not free memory",
    NULL
};


/*
 * Error handling functions
 */


/* Valid non-function to pass to libkdm functions */
static void
km_onerr_nil(int err, char *file, int line)
{
    (void) (err);
    (void) (file);
    (void) (line);
}

/* Valid non-function to pass to libkdm functions */
static void
km_onerr_print(int err, char *file, int line)
{
    fprintf(stderr, "[%s: %d] %d: %s\n", file, line, err, km_err_msgs[err]);
}

/*
 * Memory allocation/deallocation
 */

static inline void *
km_calloc_(size_t n, size_t size, void (*onerr)(int, char *, int), char *file, int line)
{
    void * ret = calloc(n, size);
    if (ret == NULL) {
        (*onerr)(1, file, line);
        return NULL;
    } else {
        return ret;
    }
}
#define km_calloc(n, sz, fn) km_calloc_(n, sz, fn, __FILE__, __LINE__)

static inline void *
km_malloc_(size_t size, void (*onerr)(int, char *, int), char *file, int line)
{
    void * ret = malloc(size);
    if (ret == NULL) {
        (*onerr)(1, file, line);
        return NULL;
    } else {
        return ret;
    }
}
#define km_malloc(sz, fn) km_malloc_(sz, fn, __FILE__, __LINE__)

static inline void *
km_realloc_(void *data, size_t size, void (*onerr)(int, char *, int), char *file, int line)
{
    void * ret = realloc(data, size);
    if (ret == NULL) {
        (*onerr)(1, file, line);
        return NULL;
    } else {
        return ret;
    }
}
#define km_realloc(ptr, sz, fn) km_realloc_(ptr, sz, fn, __FILE__, __LINE__)

#define km_free(data, onerr)        \
    STMT_BEGIN                      \
    if (data != NULL) {             \
        free(data);                 \
        data = NULL;                \
    }                               \
    STMT_END

#endif /* KDM_H */
