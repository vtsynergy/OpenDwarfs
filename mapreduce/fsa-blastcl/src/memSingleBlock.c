// memSingleBlock.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// These routines provide functionality for using a SINGLE large block of memory
// for holding individual data entries and allowing access to the portions of memory for each
// entry individually. Once an entry has been created, it CAN be moved (ie. realloced)
// however all of the entries will always be stored in a single array which can
// be used for sorting functions and alike.
// Functionality is provided to allow:
// - Getting new entries one at a time
// - Traversing all of the existing entries
// - Freeing memory used by the entire structure

#include <stdlib.h>
#include "blast.h"

// Create the memory block with an initial block size
struct memSingleBlock* memSingleBlock_initialize(size_t entrySize, int4 blockSize)
{
	struct memSingleBlock* memSingleBlock;

    // Declare memory for first block
    memSingleBlock = (struct memSingleBlock*)global_malloc(sizeof(struct memSingleBlock));

    if (memSingleBlock == NULL)
    {
		fprintf(stderr, "Error allocating memory\n");
		exit(-1);
    }

    memSingleBlock->blockSize = blockSize;
    memSingleBlock->entrySize = entrySize;
	memSingleBlock->numEntries = 0;

    // Declare memory for the block
    memSingleBlock->block = (void*)global_malloc(memSingleBlock->entrySize * memSingleBlock->blockSize);

    if (memSingleBlock->block == NULL)
    {
		fprintf(stderr, "Error allocating memory\n");
		exit(-1);
    }

    return memSingleBlock;
}

// Create the memory block with an initial block size. The struct to hold the
// memory block information is already declared and passed to this function
void memSingleBlock_initializeExisting(struct memSingleBlock* memSingleBlock,
                                       size_t entrySize, int4 blockSize)
{
    memSingleBlock->blockSize = blockSize;
    memSingleBlock->entrySize = entrySize;
	memSingleBlock->numEntries = 0;

    // Declare memory for the block
    memSingleBlock->block = (void*)global_malloc(memSingleBlock->entrySize * memSingleBlock->blockSize);

    if (memSingleBlock->block == NULL)
    {
		fprintf(stderr, "Error allocating memory\n");
		exit(-1);
    }
}

// Get an unused entry from the block
void* memSingleBlock_newEntry(struct memSingleBlock* memSingleBlock)
{
	void* newEntry;

	// Check if we need to increase the block size
	if (memSingleBlock->numEntries >= memSingleBlock->blockSize)
	{
    	// Increase the size
        memSingleBlock->blockSize *= 2;

        memSingleBlock->block = (void*)global_realloc(memSingleBlock->block,
                                 memSingleBlock->entrySize * memSingleBlock->blockSize);

        if (memSingleBlock->block == NULL)
        {
            fprintf(stderr, "Error allocating memory\n");
            exit(-1);
        }
    }

    // Use the next available slot in the latest block
    newEntry = ((char*)(memSingleBlock->block)) + memSingleBlock->numEntries * memSingleBlock->entrySize;

    memSingleBlock->numEntries++;

    return newEntry;
}

// Reset the current position to the beginning
void memSingleBlock_resetCurrent(struct memSingleBlock* memSingleBlock)
{
	memSingleBlock->currentEntry = 0;
}

// Get the current entry and advance to the next
void* memSingleBlock_getCurrent(struct memSingleBlock* memSingleBlock)
{
	void* entry;

    // If we have reached the last entry return NULL
    if (memSingleBlock->currentEntry >= memSingleBlock->numEntries)
    {
    	return NULL;
    }
	else
    {
        entry = ((char*)(memSingleBlock->block)) +
                         memSingleBlock->currentEntry * memSingleBlock->entrySize;
		memSingleBlock->currentEntry++;
        return entry;
    }
}

// Get a specific entry in the block
void* memSingleBlock_getEntry(struct memSingleBlock* memSingleBlock, int4 position)
{
	void* entry;
	entry = ((char*)(memSingleBlock->block)) + position * memSingleBlock->entrySize;
	return entry;
}

// Get the last entry in the block
void* memSingleBlock_getLastEntry(struct memSingleBlock* memSingleBlock)
{
	return ((char*)(memSingleBlock->block)) +
            (memSingleBlock->numEntries - 1) * memSingleBlock->entrySize;
}

// Free memory used by the memSingleBlock then the memSingleBlock itself
void memSingleBlock_free(struct memSingleBlock* memSingleBlock)
{
	// Free the memory block then struct itself
    free(memSingleBlock->block);
    free(memSingleBlock);
}
