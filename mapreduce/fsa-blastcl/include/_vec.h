/* _vec.h declares the implementation details of vec.[ch]
 *
 * written nml 2003-03-11
 *
 */

#ifndef PRIVATE_VEC_H
#define PRIVATE_VEC_H

#ifdef _cplusplus
extern "C" {
#endif

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

/* how much vectors grow by in vec_expand */
#define GROWTHFACTOR 2

/* byte vector structure for vbyte coding.  Note that all quantities
 * below are in bytes (not bits like the old code) */
struct vec {
    char *vector;                /* the actual bytes to hold the vector */
    unsigned long int size;      /* the capacity of vector */
    unsigned long int pos;       /* where we are in vector */
    unsigned long int len;       /* the current length of vector */
    int err;                     /* the last error to occur to the vector */
};

/* a couple of methods to initialise vectors in-place, to avoid (expensive)
 * memory management.  After being initialised by these methods, you need to
 * free() the vector component to avoid memory leaks */

int vec_init_inplace(unsigned long int initsize, struct vec *vec);

int vec_read_inplace(int srcfd, struct vec *vec, unsigned int size);

#define VEC_DELETE_INPLACE(vec) free((vec)->vector)

/* enlarge vector to size size, returns true on success and 0 on failure */
int vec_expand(struct vec* vec, unsigned int size);

#define VEC_LENGTH(vec) ((vec)->len)

#define VEC_ERR(vec) ((vec)->err)

#ifdef __MY_EXTERN_C__
}
#endif

#ifdef _cplusplus
}
#endif

#endif

