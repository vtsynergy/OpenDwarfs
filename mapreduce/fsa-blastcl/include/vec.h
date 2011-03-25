/* vec.h declares an interface to create and access bytevectors, which
 * can be dynamically read from, written to and expanded.
 *
 * If you don't know what a variable-byte integer encoding scheme is,
 * you should read Williams and Zobel, Compressing Integers for Fast
 * File Access (http://www.seg.rmit.edu.au/research/download.php?manuscript=5)
 *
 * based very closely on the old vec.[ch] written by Hugh Williams
 *
 * updated nml 2002-02-28 
 *
 */

#ifndef VEC_H
#define VEC_H

#ifdef _cplusplus
extern "C" {
#endif

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

struct vec;

/* initialise a vector and return a vector, with size initsize.  Returns NULL 
 * on failure */
struct vec* vec_init(unsigned long int initsize);

/* read a vector off disk, returns true on success and 0 on failure */
struct vec* vec_read(int srcfd, unsigned int size);

/* copy a vector (exactly, including pos but not errno), 
 * returns NULL on failure */
struct vec* vec_copy(const struct vec* vec);

/* reinitialise a vector to size initsize.  Will return a new
 * allocated vector if vec passed is NULL (think realloc).  Returns
 * NULL on failure.  If you passed an initialised vector in vec and it
 * fails then vec will be untouched, and must be deallocated as usual. */
struct vec* vec_reinit(struct vec* v, unsigned long int initsize);

/* free an initialised vector.  Returns true on success and 0 on
 * failure. */
int vec_free(struct vec* vec);

/* reset the stateful components of the vector */
void vec_resetpos(struct vec* vec);
void vec_reseterr(struct vec* vec);

/* variable byte encode n into the vector.  Returns number of bytes
 * written or 0 on failure.  Truncates the length of vector to last
 * written place. */
unsigned long int vec_putvbyte(struct vec* vec, unsigned long int n);

/* get a variable byte encoded number from the vector, and place it
 * into n.  Returns number of bytes read or 0 on failure */
unsigned long int vec_getvbyte(struct vec* vec, unsigned long int* n);

/* skip over a variable byte encoded number in the vector, returns
 * number of bytes skipped or 0 on failure */
unsigned long int vec_scanvbyte(struct vec* vec);

/* returns the length of a number as a vbyte (in bytes) */
unsigned int vec_vbyte_len(unsigned long int n);

/* return the last error that occurred to the vector (same numbers as errno).  
 * Note that hitting the end of vector will have error code ENOSPC returned */
int vec_err(const struct vec* vec);

/* returns true if the vector is positioned at the end */
int vec_eof(const struct vec* vec);

/* return the length of the vector */
unsigned long int vec_length(const struct vec* vec);

/* return the current position of the vector */
unsigned long int vec_position(const struct vec* vec);

/* write the vector, from pos to the end, to disk. Returns number of bytes
 * written on success and 0 on failure */
unsigned int vec_write(struct vec* src, int dstfd);

/* append the remaining contents of a vector (from current position to end) to
 * the end (after current position, overwriting current contents) of another.
 * returns number of bytes appended or 0 on failure. */
unsigned int vec_append(struct vec* dst, const struct vec* src);

/* compare two vectors, first by length, then by contents */
int vec_cmp(const struct vec *one, const struct vec *two);

#ifdef __MY_EXTERN_C__
}
#endif

#ifdef _cplusplus
}
#endif

#endif

