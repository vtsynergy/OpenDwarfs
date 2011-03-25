// memBlocks.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// These routines provide functionality for declaring a collection of large blocks of memory
// for holding individual data entries and allowing access to the portions of memory for each
// entry individually. Once an entry has been created, it will not be moved (ie. realloced).
// Functionality is provided to allow:
// - Getting new entries one at a time
// - Traversing all of the existing entries
// - Freeing memory used by the entire structure
// The advantage of using these routines is that they provide better locality and
// efficiency than individual malloc operations.

#include "blast.h"

// Create a new memory block
struct memBlocks* memBlocks_initialize(size_t entrySize, int4 blockSizes)
{
	struct memBlocks* memBlocks;

    // Declare memory for first block
    memBlocks = (struct memBlocks*)global_malloc(sizeof(struct memBlocks));

    memBlocks->blockSizes = blockSizes;
    memBlocks->entrySize = entrySize;
    memBlocks->numBlocks = 0;
	memBlocks->numTotalEntries = 0;
    memBlocks->maxNumBlocks = 10;

    // Declare memory for the first block
    memBlocks->lastBlock = (void*)global_malloc(memBlocks->entrySize * memBlocks->blockSizes);

    // Declare memory for pointers to blocks and add the first block
	memBlocks->blocks = (void**)global_malloc(sizeof(void*) * memBlocks->maxNumBlocks);
	memBlocks->numEntries = (int4*)global_malloc(sizeof(int4) * memBlocks->maxNumBlocks);

    memBlocks->blocks[memBlocks->numBlocks] = memBlocks->lastBlock;
    memBlocks->numEntries[memBlocks->numBlocks] = 0;
	memBlocks->numBlocks++;

    return memBlocks;
}

// Get an unused entry from the block
void* memBlocks_newEntry(struct memBlocks* memBlocks)
{
	void* newEntry;

	// Check if we need to create a new block of memory
	if (memBlocks->numEntries[memBlocks->numBlocks - 1] >= memBlocks->blockSizes)
	{
        // Declare memory for the new block
        memBlocks->lastBlock = (void*)global_malloc(memBlocks->entrySize * memBlocks->blockSizes);

        // Check if we need more memory for block pointers
        if (memBlocks->numBlocks >= memBlocks->maxNumBlocks)
        {
        	// Allocate more
            memBlocks->maxNumBlocks *= 2;
			memBlocks->blocks = (void**)global_realloc(memBlocks->blocks,
                               sizeof(void*) * memBlocks->maxNumBlocks);
            memBlocks->numEntries = (int4*)global_realloc(memBlocks->numEntries,
                               sizeof(int4) * memBlocks->maxNumBlocks);
        }

        // Store the address of this new block
        memBlocks->blocks[memBlocks->numBlocks] = memBlocks->lastBlock;

        // Reset number of entries in this block
		memBlocks->numEntries[memBlocks->numBlocks] = 0;
		memBlocks->numBlocks++;
    }

    // Use the next available slot in the latest block
    newEntry = ((char*)(memBlocks->lastBlock))
             + memBlocks->numEntries[memBlocks->numBlocks - 1] * memBlocks->entrySize;

    memBlocks->numEntries[memBlocks->numBlocks - 1]++;
    memBlocks->numTotalEntries++;
    return newEntry;
}

// Get a run of consecutive entries from the block
void* memBlocks_newEntries(struct memBlocks* memBlocks, uint4 numNewEntries)
{
	void* newEntry;

	// Check if we need to create a new block of memory
	if (memBlocks->numEntries[memBlocks->numBlocks - 1] + numNewEntries > memBlocks->blockSizes)
	{
        // Declare memory for the new block
        memBlocks->lastBlock = (void*)global_malloc(memBlocks->entrySize * memBlocks->blockSizes);

        // Check if we need more memory for block pointers
        if (memBlocks->numBlocks >= memBlocks->maxNumBlocks)
        {
        	// Allocate more
            memBlocks->maxNumBlocks *= 2;
			memBlocks->blocks = (void**)global_realloc(memBlocks->blocks,
                               sizeof(void*) * memBlocks->maxNumBlocks);
            memBlocks->numEntries = (int4*)global_realloc(memBlocks->numEntries,
                               sizeof(int4) * memBlocks->maxNumBlocks);
        }

        // Store the address of this new block
        memBlocks->blocks[memBlocks->numBlocks] = memBlocks->lastBlock;

        // Reset number of entries in this block
		memBlocks->numEntries[memBlocks->numBlocks] = 0;
		memBlocks->numBlocks++;
    }

    // Use the next available slot in the latest block
    newEntry = ((char*)(memBlocks->lastBlock))
             + memBlocks->numEntries[memBlocks->numBlocks - 1] * memBlocks->entrySize;

    memBlocks->numEntries[memBlocks->numBlocks - 1] += numNewEntries;
    memBlocks->numTotalEntries += numNewEntries;

    return newEntry;
}

// Return unused entries back to memBlocks for later allocation
void memBlocks_returnUnused(struct memBlocks* memBlocks, int4 numUnused)
{
	memBlocks->numEntries[memBlocks->numBlocks - 1] -= numUnused;
    if (memBlocks->numEntries[memBlocks->numBlocks - 1] < 0)
    {
    	fprintf(stderr, "Error: to many unused entries returned\n"); fflush(stderr);
        exit(-1);
    }
}

// Get the last entry
void* memBlocks_getLastEntry(struct memBlocks* memBlocks)
{
	void* lastEntry;

    fflush(stdout);

    if (memBlocks->numEntries[memBlocks->numBlocks - 1] < 1)
    	return NULL;

	lastEntry = ((char*)(memBlocks->lastBlock))
              + (memBlocks->numEntries[memBlocks->numBlocks - 1] - 1) * memBlocks->entrySize;

    return lastEntry;
}

// Reset the current position to the beginning
void memBlocks_resetCurrent(struct memBlocks* memBlocks)
{
	memBlocks->currentEntry = 0;
	memBlocks->currentBlock = 0;
}

// Get the current entry and advance to the next
void* memBlocks_getCurrent(struct memBlocks* memBlocks)
{
	void* entry;

    // Advance to the next block if neccessary
    if (memBlocks->currentEntry >= memBlocks->numEntries[memBlocks->currentBlock])
    {
    	memBlocks->currentBlock++;
        memBlocks->currentEntry = 0;

        // If that was the last block, return NULL
        if (memBlocks->currentBlock >= memBlocks->numBlocks)
            return NULL;
    }

    entry = ((char*)(memBlocks->blocks[memBlocks->currentBlock])) +
                     memBlocks->currentEntry * memBlocks->entrySize;
    memBlocks->currentEntry++;
    return entry;
}

// Free memory used by the memBlocks then the memBlocks itself
void memBlocks_free(struct memBlocks* memBlocks)
{
	// Free each block
	while (memBlocks->numBlocks > 0)
    {
    	memBlocks->numBlocks--;
    	free(memBlocks->blocks[memBlocks->numBlocks]);
    }

    // Free array of pointers and block sizes
    free(memBlocks->blocks);
    free(memBlocks->numEntries);

    // Free memBlocks itself
    free(memBlocks);
}
