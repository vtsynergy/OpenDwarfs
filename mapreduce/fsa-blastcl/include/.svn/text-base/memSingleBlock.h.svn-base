#ifndef _memSingleBlock_
#define _memSingleBlock_

#ifdef __cplusplus
extern "C" {
#endif

struct memSingleBlock
{
    int4 blockSize;
	size_t entrySize;
    int4 numEntries;
    void* block;
    int4 currentEntry;
};

// Create the memory block with an initial block size
struct memSingleBlock* memSingleBlock_initialize(size_t entrySize, int4 blockSize);

// Create the memory block with an initial block size. The struct to hold the
// memory block information is already declared and passed to this function
void memSingleBlock_initializeExisting(struct memSingleBlock* memSingleBlock,
                                       size_t entrySize, int4 blockSize);

// Get an unused entry from the block
extern inline void* memSingleBlock_newEntry(struct memSingleBlock* memSingleBlock);

// Reset the current position to the beginning
extern inline void memSingleBlock_resetCurrent(struct memSingleBlock* memSingleBlock);

// Get the current entry and advance to the next
extern inline void* memSingleBlock_getCurrent(struct memSingleBlock* memSingleBlock);

// Get a specific entry in the block
extern inline void* memSingleBlock_getEntry(struct memSingleBlock* memSingleBlock, int4 position);

// Get the last entry in the block
extern inline void* memSingleBlock_getLastEntry(struct memSingleBlock* memSingleBlock);

// Free memory used by the memSingleBlock then the memSingleBlock itself
void memSingleBlock_free(struct memSingleBlock* memSingleBlock);

#ifdef __cplusplus
}
#endif

#endif

