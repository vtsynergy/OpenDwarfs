// postings.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Implementation of a postings list/inverted index using a hashtable

#include "blast.h"

struct postingsList
{
	uint8 hash64bit;
	unsigned char* list;
	uint4 length;
    uint4 numEntries;
    uint4 malloced;
    uint4 lastSequenceNumber;
    struct postingsList* next;
};

uint4 postings_size;
uint4 postings_numLists;
uint4 postings_seed;
uint4 postings_hashtableSize;
struct postingsList** postings_hashtable;
unsigned char postings_vbyteEncoded[8];

// Add an entry to the postings list
void postings_addPosting(struct postingsList* postingsList, uint4 sequenceNumber, uint4 offset);
// 32-bit hashing function
uint4 postings_hash32bit(char *sequence, uint4 sequenceLength, uint4 tableSize, uint4 seed);
// 64-bit hashing function
uint8 postings_hash64bit(char *sequence, uint4 sequenceLength, uint4 seed);

// Initialize postings hashtable
void postings_initialize()
{
	uint4 count;

	postings_hashtableSize = pow(2,20);
    postings_hashtable = (struct postingsList**)global_malloc(sizeof(struct postingsList*)
                       * postings_hashtableSize);
	postings_seed = rand();
    postings_size = sizeof(struct postingsList*) * postings_hashtableSize;
    postings_numLists = 0;

    count = 0;
    while (count < postings_hashtableSize)
    {
		postings_hashtable[count] = NULL;
    	count++;
    }
}

// Add an entry to the hashtable
void postings_addEntry(unsigned char* word, uint4 wordLength, uint4 sequenceNumber, uint4 offset)
{
	struct postingsList* postingsList;
	uint4 hash32bit;
    uint8 hash64bit;

    // Calculate 32 and 64 bit hash values
	hash32bit = postings_hash32bit(word, wordLength, postings_hashtableSize, postings_seed);
	hash64bit = postings_hash64bit(word, wordLength, postings_seed);

    // Look for matching entry in hashtable
    postingsList = postings_hashtable[hash32bit];
    while (postingsList != NULL)
    {
    	if (postingsList->hash64bit == hash64bit)
        {
        	// Add new entry to postings list
			postings_addPosting(postingsList, sequenceNumber, offset);
        	return;
        }
        postingsList = postingsList->next;
    }

    // Create a new entry
    postingsList = (struct postingsList*)global_malloc(sizeof(struct postingsList));
	postingsList->hash64bit = hash64bit;
	postingsList->list = NULL;
    postingsList->length = 0;
    postingsList->malloced = 0;
    postingsList->numEntries = 0;
    postingsList->lastSequenceNumber = 0;
	postings_size += sizeof(struct postingsList);
    postings_numLists++;

    // Add to start of list
	postingsList->next = postings_hashtable[hash32bit];
    postings_hashtable[hash32bit] = postingsList;

    // Add new entry to postings list
	postings_addPosting(postingsList, sequenceNumber, offset);
}

// Compare the length of two entries in the hashtable
int4 postings_compareLists(const void* list1, const void* list2)
{
	const struct finalList *e1, *e2;

	e1 = (struct finalList*)list1;
	e2 = (struct finalList*)list2;

	if (e1->numEntries > e2->numEntries)
	{
		return -1;
	}
	else if (e1->numEntries < e2->numEntries)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Return sorts lists and free postings hash structure in the process
struct finalList* postings_getSortedLists()
{
	struct finalList* finalLists;
    struct postingsList* postingsList, *previousPostingsList;
    uint4 listCount = 0, totalSize = 0;
	uint4 hash32bit;

    finalLists = (struct finalList*)global_malloc(sizeof(struct finalList) * postings_numLists);

    // For each list
	hash32bit = 0;
    while (hash32bit < postings_hashtableSize)
    {
        postingsList = postings_hashtable[hash32bit];
        while (postingsList != NULL)
        {
        	// Add to new list of final lists
			finalLists[listCount].numEntries = postingsList->numEntries;
			finalLists[listCount].list = postingsList->list;

			totalSize += postingsList->length;

			listCount++;
            previousPostingsList = postingsList;
            postingsList = postingsList->next;

            free(previousPostingsList);
        }

    	hash32bit++;
    }

    free(postings_hashtable);

    // Sort the lists in order of number of entries
	qsort(finalLists, postings_numLists, sizeof(struct finalList), postings_compareLists);

    // Remove lists with less than 2 entries
    while (postings_numLists > 0 && finalLists[postings_numLists - 1].numEntries < 2)
    {
		postings_numLists--;
        free(finalLists[postings_numLists].list);
    }

//	printf("Total postings size=%d\n", totalSize);

    return finalLists;
}

// Decode a list and free at the same time
uint4 postings_decodeList(struct finalList finalList, struct wordLocation* wordLocations)
{
	unsigned char* list;
	uint4 numEntries, sGap, offset, sequenceNumber = 0, entry = 0;

    list = finalList.list;
    numEntries = finalList.numEntries;

    while (entry < numEntries)
    {
        vbyte_getVbyte(list, &sGap);
        sequenceNumber += sGap;
        vbyte_getVbyte(list, &offset);

		wordLocations[entry].sequenceNumber = sequenceNumber;
		wordLocations[entry].offset = offset;

//        printf("%d,%d (%d/%d entries)\n", sequenceNumber, offset, entry, numEntries); fflush(stdout);

    	entry++;
    }

    free(finalList.list);

    return numEntries;
}

// Add an entry to the postings list
void postings_addPosting(struct postingsList* postingsList, uint4 sequenceNumber, uint4 offset)
{
	uint4 sGap, encodedLength;
    unsigned char* vbyteEncoded;

    // Ignore multiple hits for same sequence
    if (sequenceNumber == postingsList->lastSequenceNumber && sequenceNumber != 0)
    	return;

    // Vbyte encode sequence number, offset information
    sGap = sequenceNumber - postingsList->lastSequenceNumber;
    vbyteEncoded = postings_vbyteEncoded;
    vbyte_putVbyte(vbyteEncoded, sGap);
    vbyte_putVbyte(vbyteEncoded, offset);
    encodedLength = vbyteEncoded - postings_vbyteEncoded;

    // If first entry in list
	if (postingsList->list == NULL)
    {
    	// Initialize list
    	postingsList->malloced = encodedLength;
    	postingsList->list = global_malloc(postingsList->malloced);
        memcpy(postingsList->list, postings_vbyteEncoded, encodedLength);
        postingsList->length = encodedLength;
		postings_size += postingsList->malloced;
    }
    else
    {
    	// If list is not large enough
    	if (encodedLength + postingsList->length > postingsList->malloced)
        {
        	// Increase list size
            postings_size -= postingsList->malloced;
			postingsList->malloced = postingsList->malloced * 2 + encodedLength;
			postingsList->list = global_realloc(postingsList->list, postingsList->malloced);
            postings_size += postingsList->malloced;
        }

        // Append data
		memcpy(postingsList->list + postingsList->length, postings_vbyteEncoded, encodedLength);
        postingsList->length += encodedLength;
    }

    postingsList->lastSequenceNumber = sequenceNumber;
    postingsList->numEntries++;
}

// Print postings lists
void postings_print()
{
	printf("%d lists. Postings size=%.2f Mb\n", postings_numLists, (float)postings_size / 1024.0 / 1024.0);
}

// 32-bit hashing function
uint4 postings_hash32bit(char *sequence, uint4 sequenceLength, uint4 tableSize, uint4 seed)
{
   unsigned int c;
   unsigned long hashval = seed;

   for (c = 0; c < sequenceLength; c++)
   {
         hashval ^= ((hashval << 7) + sequence[c] + (hashval >> 2));
   }

   return hashval % tableSize;
}

// 64-bit hashing function
uint8 postings_hash64bit(char *sequence, uint4 sequenceLength, uint4 seed)
{
   int c;
   uint8 hashval = seed;

   for (c = 0; c < sequenceLength; c++)
         hashval ^= ((hashval << 5) + sequence[c] + (hashval >> 2) + (hashval << 11));

   return hashval;
}

