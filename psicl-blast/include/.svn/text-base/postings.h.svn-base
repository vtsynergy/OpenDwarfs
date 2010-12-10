#ifndef _postings_
#define _postings_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

extern uint4 postings_numLists;

struct wordLocation
{
	uint4 sequenceNumber;
    uint4 offset;
};

struct finalList
{
	uint4 numEntries;
	unsigned char* list;
};

// Initialize postings hashtable
void postings_initialize();

// Add an entry to the hashtable
void postings_addEntry(unsigned char* word, uint4 wordLength, uint4 sequenceNumber, uint4 offset);

// Print postings lists
void postings_print();

// Return sorts lists and free postings hash structure in the process
struct finalList* postings_getSortedLists();

// Decode a list and free at the same time
uint4 postings_decodeList(struct finalList finalList, struct wordLocation* wordLocations);

#ifdef __MY_EXTERN_C__
}
#endif

#endif

