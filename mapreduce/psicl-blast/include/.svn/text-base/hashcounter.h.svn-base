/** @file
 *
 * The hashcounter is an efficient hash-based phrase occurrence counter
 * specifically designed for the SPEX (Shared Phrase EXtraction) algorithm.
 * The only requirement from the SPEX algorithm of this data structure is
 * that it maintain a record or whether a given phrase occurred in the
 * collection or not. Because false positives are acceptable but false
 * negatives are not, a basic hash structure of an appropriate size is a
 * good choice.
 *
 * Because the only information that is required is whether a phrase
 * appears multiple times (ie. the exact count is irrelevant), two bits per
 * entry are sufficient. This allows for a large hashtable to be stored
 * quite compactly, minimising the probability of collisions.
 *
 * $Id: hashcounter.h,v 1.12.4.1 2005/05/05 07:24:07 ybernste Exp $
 * Written by Yaniv Bernstein 2004
*/

#ifndef HASHCOUNTER_H_
#define HASHCOUNTER_H_

#include <stdio.h>

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

/** ADT for the hashcounter. Contains storage and metadata */
typedef struct
{
   char *hashtable; /**< Block of memory for the actual hashtable */
   unsigned long size; /**< The size of the hashtable */
   unsigned long byte_size; /**< The size of the hashtable */
   unsigned long indmask; 
   unsigned long fieldmask; 
   unsigned long hashseed; /**< The seed for the hashing process */
   /** The number of insertions made into the hashcounter since it was last
    * reset. */
   unsigned long insertions;
} Hashcounter;

/** Initialise a new hashcounter of size size. Note that size is the size in
 * bytes -- the hashtable will in fact contain 4 * bytes fields */
Hashcounter *hashcounter_new(unsigned long size, unsigned long seed);

/** Increment the count for the given phrase within the hashcounter. The
 * phrase is given as an array of strings, each containing a single word.
 * 
 * The parameter 'words' indicates the number of words in the phrase, while
 * 'array_size' indicates the size of the array (which is sometimes bigger
 * when hashing subphrases). The 'offset' parameter indicates the position
 * in the array in which the phrase starts -- the array is circular, so
 * when the end of the array is reached, the phrase is assumed to continue
 * at the beginning of the array. */
void hashcounter_insert(Hashcounter *counter, uint4 hashValue);

/** Mainly for internal use: return a hash value for a given phrase. See the
 * comment for 'hashcounter_new' for notes on the meanings of the various
 * parameters. */
//unsigned long hashcounter_hash(char *chunk, int chunk_length);


/** Returns a nonzero value value if the phrase in question appears more
 * than once in the hashcounter. See the comment for 'hashcounter_new' for
 * notes on the meanings of the various parameters. */
/*int hashcounter_multiple(Hashcounter *counter, char **phrase, int word_count, int array_size, int offset);*/
int hashcounter_multiple(Hashcounter *counter, uint4 hashValue);

int hashcounter_single(Hashcounter *counter, uint4 hashValue);

/** returns the number of phrases that have been inserted into this
 * hashcounter */
unsigned long hashcounter_insertions(Hashcounter *counter);

/** Write the hashtable out to stream 'stream' */
void hashcounter_write(Hashcounter *counter, FILE *stream);

/** Print out a count of fields with different values */
void hashcounter_count(Hashcounter *counter, FILE *stream);

/** Clear the hashcounter */
void hashcounter_reset(Hashcounter *counter);

/** Clean up the hascounter when you're done. */
void hashcounter_free(Hashcounter *counter);

// Hashing function
#define hashcounter_hash(counter, chunk, chunk_length, hashValue, position) \
if (1) \
{ \
    hashValue = 5381; \
    for (position = 0; position < chunk_length; position++) \
    { \
          hashValue ^= ((hashValue << 7) + chunk[position] + (hashValue >> 2)); \
    } \
}

#ifdef __MY_EXTERN_C__
}
#endif

#endif
