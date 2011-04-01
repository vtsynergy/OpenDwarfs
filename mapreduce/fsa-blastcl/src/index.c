// index.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for creating an in-memory inverted index

#include "blast.h"

struct indexCoordinate** index_sequenceCoordinates;
struct indexCoordinate* index_coordinates;
uint4 index_numCoordinates;

uint4* index_sequencePositions;
uint4* index_descriptionLocations;

struct wordList* index_words;
uint4 index_numWords;
unsigned char index_vbyteEncoded[4];
uint4 index_currentCoordinate;

uint4* index_loadedWords;
unsigned char* index_offsets;
uint4 index_subjectNumber;

uint4 index_wordSize = 12, index_intervalSize = 9;
#define index_charSize 2

// Compare the two queryWords' codewords
int4 alignments_compareCodeword(const void* queryWord1, const void* queryWord2);
// Compare the two queryWords' queryPositions
int4 alignments_compareQueryPosition(const void* queryWord1, const void* queryWord2);
void index_addWord(uint4 codeword, uint4 subjectNumber, uint4 offset);
uint4 index_generateCodeword(unsigned char* word, uint4 wordSize);

// **** INDEX CREATION CODE ****

// Initialize the creation of a new index structure
void index_initializeBuild(uint4 fromCodeword, uint4 toCodeword)
{
	uint4 codeword;

//	index_numWords = pow(4, index_wordSize);
    index_words = (struct wordList*)global_malloc(sizeof(struct wordList) * (toCodeword - fromCodeword));
	index_words -= fromCodeword;

    // For each word
    codeword = fromCodeword;
    while (codeword < toCodeword)
    {
    	// Initialize list of occurrences
    	index_words[codeword].offsets = NULL;
		index_words[codeword].length = 0;
        index_words[codeword].allocated = 0;
        index_words[codeword].lastOffset = 0;
        index_words[codeword].lastSequenceNumber = 0;

    	codeword++;
    }

    index_subjectNumber = 0;
}

// Finish building inverted lists for range of codewords
void index_finishBuild(uint4 fromCodeword, uint4 toCodeword)
{
	uint4 codeword;

    // For each word
    codeword = fromCodeword;
    while (codeword < toCodeword)
    {
    	// Free list
    	free(index_words[codeword].offsets);
    	codeword++;
    }

    // Free the lists
	index_words += fromCodeword;
    free(index_words);
}

// Index every Nth word of length index_wordSize in the subject
void index_addSubject(unsigned char* subject, uint4 subjectLength, uint4 fromCodeword, uint4 toCodeword)
{
	uint4 codeword, wordOffset;

    // Add subjectNumber and offset to lists for each word
	wordOffset = 0;
    while (wordOffset + index_wordSize - 1 < subjectLength)
    {
    	// For this word
    	codeword = index_generateCodeword(subject + wordOffset, index_wordSize);

        // If it is in the range of codewords we are considering
        if (codeword >= fromCodeword && codeword < toCodeword)
        {
            // Record subject number and offset of word in index
            index_addWord(codeword, index_subjectNumber, wordOffset / index_intervalSize);// + (index_wordSize - 4));
		}

        // Index every Nth word
        wordOffset += index_intervalSize;
    }

    index_subjectNumber++;
}

// Add a given word position to the index
void index_addWord(uint4 codeword, uint4 subjectNumber, uint4 offset)
{
	struct wordList* wordList;
    unsigned char* vbyteEncoded;
    uint4 sequenceGap, offsetGap, vbyteSize;

    wordList = index_words + codeword;

    // Encoded the difference between this sequence number and the last
    sequenceGap = subjectNumber - wordList->lastSequenceNumber;
    vbyteEncoded = index_vbyteEncoded;
	vbyte_putVbyte(vbyteEncoded, sequenceGap);

    // If last occurrence was in the same sequence, store difference in offet
    if (sequenceGap == 0)
    	offsetGap = offset - wordList->lastOffset;
	else
	    // Else store absolute offset
		offsetGap = offset;
    vbyte_putVbyte(vbyteEncoded, offsetGap);

    vbyteSize = vbyteEncoded - index_vbyteEncoded;

    // If we need more space in the list of offsets
    if (wordList->length + vbyteSize >= wordList->allocated)
    {
		wordList->allocated += vbyteSize + 10;

        wordList->offsets = (unsigned char*)global_realloc(wordList->offsets,
                            sizeof(char) * wordList->allocated);
    }

	// Record encoded offset
    memcpy(wordList->offsets + wordList->length, index_vbyteEncoded, vbyteSize);

    wordList->lastOffset = offset;
    wordList->lastSequenceNumber = subjectNumber;
    wordList->length += vbyteSize;
}

// Return number of hit offsets for given codeword
uint4 index_numWordOffsets(uint4 codeword)
{
	return index_words[codeword].length;
}

// Return the list of hit offsets for given codeword
unsigned char* index_wordOffsets(uint4 codeword)
{
	return index_words[codeword].offsets;
}



// **** QUERY PROCESS CODE ****

// Given a query sequence uses an inverted index of the collection to identify the
// sequence number and offset of all hits between the query and the collection
void index_processQuery(unsigned char* startIndex, struct PSSMatrix PSSMatrix,
                        uint4 numSequences)
{
	uint4 queryPosition, codeword = 0, queryPosition4;
    unsigned char* offsets, *endOffsets;
    uint4 offsetGap, offset, sequenceGap, sequenceNumber;
    struct indexCoordinate* coordinate;
	struct memBlocks* unsortedCoordinates;
    uint4 *numSubjectHits, numQueryPositions, queryWordCount, numOffsets;
    uint4 time, wordPosition, containsWildcard;
	struct queryWord* queryWords;

    // Read word and interval size from start of index
	vbyte_getVbyte(startIndex, &index_wordSize);
	vbyte_getVbyte(startIndex, &index_intervalSize);

	index_numWords = pow(4, index_wordSize);
    index_sequencePositions = (uint4*)startIndex;
    index_descriptionLocations = index_sequencePositions + numSequences;
	index_loadedWords = index_descriptionLocations + numSequences;
	index_offsets = (unsigned char*)(index_loadedWords + index_numWords + 1);

    time = clock();
    unsortedCoordinates = memBlocks_initialize(sizeof(struct indexCoordinate), numSequences);

    // Declare and initialize array for count number of hits for each sequence
    numSubjectHits = (uint*)global_malloc(sizeof(uint4) * numSequences);
	sequenceNumber = 0;
    while (sequenceNumber < numSequences)
    {
    	numSubjectHits[sequenceNumber] = 0;
    	sequenceNumber++;
    }

    // Memory to hold offsets string for each query word
    numQueryPositions = PSSMatrix.length - index_wordSize + 1;
	queryWords = (struct queryWord*)global_malloc(sizeof(struct queryWord) * numQueryPositions);

    // For each word in the query
    queryPosition = 0;
    while (queryPosition < numQueryPositions)
    {
    	// Check if the word contains a wildcard
        containsWildcard = 0; wordPosition = 0;
        while (wordPosition < index_wordSize)
        {
            if (PSSMatrix.queryCodes[queryPosition + wordPosition] >= encoding_numRegularLetters)
                containsWildcard = 1;

            wordPosition++;
        }

        // Don't include words that cross the strand boundry or contain wildcards
        if (!containsWildcard && !(queryPosition < PSSMatrix.strandLength &&
              queryPosition >= PSSMatrix.strandLength - index_wordSize + 1))
		{
//            printf("--Query position=%d\n", queryPosition);

            // Get the codeword
            codeword = index_generateCodeword(PSSMatrix.bestMatchCodes + queryPosition, index_wordSize);

            // Get wordlist for that codeword
            offsets = index_offsets + index_loadedWords[codeword];
            endOffsets = index_offsets + index_loadedWords[codeword + 1];

            queryWords[queryPosition].offsets = offsets;
			queryWords[queryPosition].endOffsets = endOffsets;
			queryWords[queryPosition].queryPosition = queryPosition;
            queryWords[queryPosition].codeword = codeword;

//            printf("codeword=%d start=%d end=%d numHits=%d\n", codeword, index_loadedWords[codeword],
//                   index_loadedWords[codeword + 1], endOffsets - offsets);
		}
        else
        {
            queryWords[queryPosition].offsets = NULL;
			queryWords[queryPosition].endOffsets = NULL;
			queryWords[queryPosition].queryPosition = queryPosition;
            queryWords[queryPosition].codeword = codeword;
        }

//        printf("\n");
    	queryPosition++;
    }

    // Sort the query words by codeword
	qsort(queryWords, numQueryPositions, sizeof(struct queryWord), alignments_compareCodeword);

    // For each query word
    queryWordCount = 0;
    while (queryWordCount < numQueryPositions)
    {
    	// Ignoring those that cross the strand boundry
		if (queryWords[queryWordCount].offsets != NULL)
        {
        	// Make in-memory copy of list of offsets
            numOffsets = queryWords[queryWordCount].endOffsets - queryWords[queryWordCount].offsets;
			offsets = (char*)global_malloc(sizeof(char) * numOffsets);

            memcpy(offsets, queryWords[queryWordCount].offsets, numOffsets);
			queryWords[queryWordCount].offsets = offsets;
            queryWords[queryWordCount].endOffsets = offsets + numOffsets;
		}

        queryWordCount++;
    }

    // Sort the query words by query position
	qsort(queryWords, numQueryPositions, sizeof(struct queryWord), alignments_compareQueryPosition);

    queryPosition = 0;
    while (queryPosition < numQueryPositions)
    {
    	// Ignoring those that cross the strand boundry
		if (queryWords[queryPosition].offsets != NULL)
        {
        	offsets = queryWords[queryPosition].offsets;
            endOffsets = queryWords[queryPosition].endOffsets;
            offset = 0;
            sequenceNumber = 0;
        	queryPosition4 = queryPosition + (index_wordSize - 4);

            // Traverse the offsets
            while (offsets < endOffsets)
            {
                vbyte_getVbyte(offsets, (&sequenceGap));
                vbyte_getVbyte(offsets, (&offsetGap));

//                printf("[%d,%d]\n", sequenceGap, offsetGap);

                if (sequenceGap > 0)
                {
                	offset = offsetGap;
                    sequenceNumber += sequenceGap;
                }
                else
                {
                	offset += offsetGap;
                }
    //            printf(" %u", offset);

                // Add query/database coordinate of match to relevant bucket
//                printf("Sequence number=%d\n", sequenceNumber);
                coordinate = (struct indexCoordinate*)memBlocks_newEntry(unsortedCoordinates);
                coordinate->queryOffset = queryPosition4;
                coordinate->subjectOffset = offset * index_intervalSize + (index_wordSize - 4);
                coordinate->subjectNumber = sequenceNumber;

                numSubjectHits[sequenceNumber]++;
//                printf("[%d,%d]\n", queryPosition, offset);

                blast_numHits++;
            }

            free(queryWords[queryPosition].offsets);
		}

        queryPosition++;
	}


    printf("Time to process query=%f\n", (float)(clock() - time) / CLOCKS_PER_SEC);
    time = clock();

    // Make memory for sorted list
    index_numCoordinates = unsortedCoordinates->numTotalEntries;
	index_coordinates = (struct indexCoordinate*)global_malloc(
                         sizeof(struct indexCoordinate) * index_numCoordinates);
	index_sequenceCoordinates = (struct indexCoordinate**)global_malloc(
                                 sizeof(struct indexCoordinate*) * numSequences);

    // For each sequence
	coordinate = index_coordinates;
    sequenceNumber = 0;
    while (sequenceNumber < numSequences)
    {
    	// If it has hits
    	if (numSubjectHits[sequenceNumber] != 0)
        {
        	// Point to location in sorted list of coordinates
			index_sequenceCoordinates[sequenceNumber] = coordinate;
            coordinate += numSubjectHits[sequenceNumber];

            numSubjectHits[sequenceNumber] = 0;
        }
    	sequenceNumber++;
    }

    // Move through list of unsorted coordinates
    memBlocks_resetCurrent(unsortedCoordinates);
    while ((coordinate = memBlocks_getCurrent(unsortedCoordinates)) != NULL)
    {
    	sequenceNumber = coordinate->subjectNumber;
//    	printf("%d,%d=[%d]\n", index_sequenceCoordinates[sequenceNumber], numSubjectHits[sequenceNumber], sequenceNumber);
    	// Place into sorted list
		index_sequenceCoordinates[sequenceNumber][numSubjectHits[sequenceNumber]] = *coordinate;
		numSubjectHits[sequenceNumber]++;
    }

    memBlocks_free(unsortedCoordinates);

/*    // Print sorted coordinates
	coordinate = index_coordinates;
    while (coordinate < index_coordinates + index_numCoordinates)
    {
    	printf("[%d]", coordinate);
    	printf("Subject %d Offset %d,%d\n", coordinate->subjectNumber, coordinate->queryOffset,
                                            coordinate->subjectOffset);
    	coordinate++;
    }*/

    printf("Time to sort buckets=%f\n", (float)(clock() - time) / CLOCKS_PER_SEC);
}

// Get the first coordinate in buckets
struct indexCoordinate* index_getFirstCoordinate()
{
    // Reset counters
	index_currentCoordinate = 0;

    // Get coordinate
    return index_getNextCoordinate();
}

// Get the next coordinate in available buckets
struct indexCoordinate* index_getNextCoordinate()
{
	struct indexCoordinate* coordinate;

    if (index_currentCoordinate >= index_numCoordinates)
    	return NULL;

    // Get current coordinate and return it
	coordinate = index_coordinates + index_currentCoordinate;
	index_currentCoordinate++;

    return coordinate;
}

// Compare the two queryWords' queryPositions
int4 alignments_compareQueryPosition(const void* queryWord1, const void* queryWord2)
{
	const struct queryWord *q1, *q2;

	q1 = (struct queryWord*)queryWord1;
	q2 = (struct queryWord*)queryWord2;

	if (q1->queryPosition > q2->queryPosition)
	{
		return 1;
	}
	else if (q1->queryPosition < q2->queryPosition)
	{
		return -1;
	}
	else
	{
    	return 0;
	}
}

// Compare the two queryWords' codewords
int4 alignments_compareCodeword(const void* queryWord1, const void* queryWord2)
{
	const struct queryWord *q1, *q2;

	q1 = (struct queryWord*)queryWord1;
	q2 = (struct queryWord*)queryWord2;

	if (q1->codeword > q2->codeword)
	{
		return 1;
	}
	else if (q1->codeword < q2->codeword)
	{
		return -1;
	}
	else
	{
    	return 0;
	}
}


// ***** SHARED CODE *****

// Generate a codeword from a given word
uint4 index_generateCodeword(unsigned char* word, uint4 wordSize)
{
	uint4 codeword = 0;

	while (wordSize > 0)
    {
		codeword = (codeword << index_charSize) | *word;
    	word++; wordSize--;
    }

    return codeword;
}

// Print the contents of the index
void index_print()
{
	uint4 codeword = 0;
	struct wordList* wordList;
    unsigned char* offsets, *endOffsets;
    uint4 offsetGap, offset, numOffsets;
	uint4 totalSize = 0;

    while (codeword < index_numWords)
    {
    	numOffsets = 0; offset = 0;

        wordList = index_words + codeword;

		offsets = wordList->offsets;
        endOffsets = offsets + wordList->length;

        totalSize += wordList->length;

        if (offsets < endOffsets)
        	printf("\nCodeword=%u:", codeword);

        while (offsets < endOffsets)
        {
			vbyte_getVbyte(offsets, (&offsetGap));
            offset += offsetGap;
            printf(" %u", offset);
            numOffsets++;
        }

//    	printf("[%d/%d = %f]\n", wordList->length, numOffsets, (float)(wordList->length) / (float)numOffsets);

        codeword++;
    }

    printf("\nTotal table size=%d bytes\n", totalSize);
}

