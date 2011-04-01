/* vbyte.h declares functions to read and write variable-byte encoded
 * integers to files
 *
 * written Hugh Williams
 *
 * updated nml 2003-03-07
 *
 */

#ifndef VBYTE_H
#define VBYTE_H

#ifdef _cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

/* read a number from a file into n, returning the number of bytes
 * read or 0 on failure */
unsigned int vbyte_read(FILE* fp, unsigned long int* n);

/* write a number n to file, returning number of bytes written or 0 on
 * failure */
unsigned int vbyte_write(FILE* fp, unsigned long int n);

/* read over a number from file, returns number of bytes scanned or 0
 * on failure */
unsigned int vbyte_scan(FILE* fp);

/* return the on disk length of the number n */
unsigned int vbyte_len(unsigned long int n);

#ifdef __MY_EXTERN_C__
}
#endif

#ifdef _cplusplus
}
#endif

#endif

