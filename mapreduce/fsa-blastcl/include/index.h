#ifndef _index_
#define _index_

#ifdef __cplusplus
extern "C" {
#endif

extern uint4 index_wordSize, index_intervalSize;

// **** INDEX CREATION CODE ****

// Initialize the creation of a new index structure
void index_initializeBuild(uint4 fromCodeword, uint4 toCodeword);
// Index every Nth word of length index_wordSize in the subject
void index_addSubject(unsigned char* subject, uint4 subjectLength, uint4 fromCodeword, uint4 toCodeword);
// Finish building inverted lists for range of codewords
void index_finishBuild(uint4 fromCodeword, uint4 toCodeword);

// **** QUERY PROCESS CODE ****

// Return an array of locations in the index file of the offsets for each word
uint4* index_wordOffsetPositions();
// Return number of hit offsets for given codeword
uint4 index_numWordOffsets(uint4 codeword);
// Return the list of hit offsets for given codeword
unsigned char* index_wordOffsets(uint4 codeword);
// Given a query sequence and inverted index of the collection, identifies the
// sequence number and offset of all hits between the query and the collection
void index_processQuery(unsigned char* startIndex, struct PSSMatrix PSSMatrix, uint4 collectionSize);
// Get the first coordinate in buckets
struct indexCoordinate* index_getFirstCoordinate();
// Get the next coordinate in available buckets
struct indexCoordinate* index_getNextCoordinate();

// **** SHARED FUNCTIONS ****

// Generate a codeword from a given word
uint4 index_generateCodeword(unsigned char* word, uint4 wordSize);
// Print the contents of the index
void index_print();

extern struct indexCoordinate** index_sequenceCoordinates;
extern struct indexCoordinate* index_coordinates;
extern uint4 index_numCoordinates;

extern uint4* index_sequencePositions;
extern uint4* index_descriptionLocations;

struct indexCoordinate
{
	uint4 queryOffset;
    uint4 subjectOffset;
    uint4 subjectNumber;
};

struct wordList
{
	unsigned char* offsets;
    uint4 length;
    uint4 allocated;
    uint4 lastOffset;
    uint4 lastSequenceNumber;
};

struct queryWord
{
	uint4 codeword;
    uint4 queryPosition;
    char* offsets, *endOffsets;
};

#ifdef __cplusplus
}
#endif

#endif

