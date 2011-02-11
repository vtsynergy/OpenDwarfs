#ifndef _memBlocks_
#define _memBlocks_

#ifdef __cplusplus
extern "C" {
#endif

struct memBlocks
{
    int4 blockSizes;	// Maximum size of each block
	size_t entrySize;	// Size of each entry
    int4 numBlocks;		// Number of blocks
    int4* numEntries;	// Number of entries for each of the blocks
    int4 maxNumBlocks;	// Maximum number of blocks before realloc
    void** blocks;		// The blocks
    void* lastBlock;	// The last block
    int4 currentEntry;	// The current entry number
    int4 currentBlock;	// The current block number
    int4 numTotalEntries;  // Total number of entries
};

// Create a new memory block
struct memBlocks* memBlocks_initialize(size_t entrySize, int4 blockSizes);

// Get an unused entry from the block
void* memBlocks_newEntry(struct memBlocks* memBlocks);

// Get a run of consecutive entries from the block
void* memBlocks_newEntries(struct memBlocks* memBlocks, uint4 numNewEntries);

// Return unused entries back to memBlocks for later allocation
void memBlocks_returnUnused(struct memBlocks* memBlocks, int4 numUnused);

// Get the last entry
void* memBlocks_getLastEntry(struct memBlocks* memBlocks);

// Reset the current position to the beginning
extern inline void memBlocks_resetCurrent(struct memBlocks* memBlocks);

// Get the current entry and advance to the next
extern inline void* memBlocks_getCurrent(struct memBlocks* memBlocks);

// Free memory used by the memBlocks then the memBlocks itself
extern inline void memBlocks_free(struct memBlocks* memBlocks);
#ifdef __cplusplus
}
#endif

#endif

