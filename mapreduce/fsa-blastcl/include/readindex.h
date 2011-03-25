/** @file
 *
 * Functions designed for sequentially extracting postings lists from lucy
 * index files, primarily designed for use with the spex indices from the 
 * deco algorithm. 
 *
 * Yaniv Bernstein 2004
 */

#ifndef READINDEX_H
#define READINDEX_H

#include <stdio.h>

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

/** Data structure for storing buffered file pointers, makes it convenient
 * when extracting vectors */
struct index_scanner
{
   FILE *idx_file; /**< A pointer to the idx file for this index */
   FILE *offsets; /**< A pointer to the file containing long lists and their offsets */
};

struct lo_info
{
   unsigned long size;
   long offset;
};

/** Stores relevant information from an extracted phrase */
struct listinfo
{
   /** Pointer to the current phrase */
   char *phrase;
   /** The number of documents the phrase occurs in */
   unsigned long doc_count;
   unsigned long size;
   /** Array of the document numbers the phrase occurs in */
   unsigned long *doc_numbers;
   /** Array of the frequency of the phrase's occurrence in each document */
   unsigned long *phrase_frequency;
   /** Array of offsets for each phrase, for each document that the phrase occurs in */
   unsigned long **phrase_offsets;
};


/** Returns a data structure that contains open buffered file pointers to
 * the index file and all the vector files in preparation for sequential
 * reading of postings vectors. */
struct index_scanner *open_index(char *index_prefix);

/** Return a phrase_info data structure for the next phrase in the index */
struct listinfo *get_next_list(struct index_scanner *iscn, struct listinfo *info);

struct listinfo *get_list_at(struct index_scanner *iscn, struct listinfo *info, long offset);

struct listinfo *get_next_biglist(struct index_scanner *iscn, struct listinfo *info);

struct lo_info *get_next_lo_info(struct index_scanner *iscn, struct lo_info *lo_info);

/** Deallocates the memory used by the data structure */
void info_free(struct listinfo *info);

/** Reset all the file pointers back to the beginning of the index */
struct index_scanner *reset_index(struct index_scanner *scn);

/** Close all the file pointers to the index */
struct index_scanner *close_index(struct index_scanner *scn);

#ifdef __MY_EXTERN_C__
}
#endif

#endif
