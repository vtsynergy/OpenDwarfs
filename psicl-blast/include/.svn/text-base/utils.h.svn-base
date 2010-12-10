/**
 * @file
 *
 * $Id: utils.h,v 1.10.2.2 2005/03/15 01:21:10 ybernste Exp $
 * Written by Yaniv Bernstein 2004
 */
#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdint.h>

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

/** A wrapper around malloc that spews if NULL is returned */
void *check_malloc(size_t size);

/** A wrapper around realloc that spews if a NULL pointer is returned */
void *check_realloc(void *ptr, size_t size);

/** A wrapper around calloc that spews if a NULL pointer is returned */
void *check_calloc(size_t elements, size_t size);

/** Takes an array of words and returns a single string composed of all
 * those words appended to each other, separated by spaces. */
char *word_flatten(char **words, int start_word, int n);

unsigned long data_hash(void *pack, unsigned long ecount, unsigned long esize, unsigned long elength, unsigned long offset, unsigned long hashsize);

uint64_t data_hash64(void *pack, unsigned long datalen);

/** Hash a string using the Karp-Rubin string hashing algorithm */
unsigned long long string_hash(char *string, unsigned long long size);

/** Print the current date and time in a nicely formatted way */
void print_time(FILE *stream);

/** Print the usage information for this program */
void usage(char *pname);

void casefold(char *string);

int fcomp(const void *f1, const void *f2);

#ifdef __MY_EXTERN_C__
}
#endif

#endif
