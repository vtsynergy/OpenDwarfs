// wordLookupDFA.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Deterministic Finite Automata for performing Wilbur-Lipman style fast lookup
// search with neighbours (non-exact word matches)

#include <stdlib.h>
#include <math.h>

#include "blast.h"

uint2* wordLookupDFA_additionalQueryPositions = NULL;
int4 wordLookupDFA_numAdditionalQueryPositions;
struct group *wordLookupDFA_groups = NULL;
//Shucai
uint4 wordLookupDFA_numWords = 0;
uint4 wordLookupDFA_numGroups = 0;
uint4 wordLookupDFA_numExtPositions = 0;
struct groupFP *wordLookupDFA_groupsFP = NULL;

int4 wordLookupDFA_numCodes, wordLookupDFA_wordLength, wordLookupDFA_numBlocks;

struct initialWord
{
	int2 numQueryPositions;
	int2 allocQueryPositions;
    int2* queryPositions;
};

struct aaFrequencyGroup
{
	struct memSingleBlock* groups;
    float frequency;
};

struct neighbour
{
	int4 codeword;
    int4 score;
    int4 position;
};

unsigned char *wordLookupDFA = NULL;
unsigned char **wordLookupDFA_blockAddresses = NULL;
//Shucai
//unsigned int *wordLookupDFA_blockAddressesFP;

void wordLookupDFA_printCodes(int4 codeword, int4 wordLength);
struct aaFrequencyGroup* wordLookupDFA_calcFrequencyGroups(struct memSingleBlock** groupsHash,
                                     struct PSSMatrix PSSMatrix, int4 numCodes, int4 wordLength);
inline int4 wordLookupDFA_getCodeword(unsigned char* codes, int4 wordLength);
inline unsigned char* wordLookupDFA_getCodes(int4 codeword, int4 wordLength);
void wordLookupDFA_getNeighbours(struct PSSMatrix PSSMatrix, int4 queryPosition,
                                 int4* numNeighbours, struct neighbour* neighbours);
void wordLookupDFA_transformationForGPU(int4 numCodes, int4 wordLength);

// Build the word-lookup structure                                 
void wordLookupDFA_build(struct PSSMatrix PSSMatrix, int4 numCodes, int4 wordLength)
{
	struct initialWord* initialWords, *initialWord;
    int4 numExternalPositions = 0;
    unsigned char *currentBlock, *currentPosition;
    int2 *startExternalPositions, *endExternalPositions;
    int2 *externalPositions[numCodes];
	int4 count, count2, groupCount, groupCodeword, match;
	struct aaFrequencyGroup* aaFrequencyGroups;
    int4 codeword, codeword2, queryPosition, *groupList1, *groupList2 = NULL;
	unsigned char code;
	int4 queryPositionCount;
    struct neighbour* neighbours = NULL;
    int4 numNeighbours;
    int4 numWords, numGroups;
	struct memSingleBlock** groupsHash, **groupLists, *groupList;
	struct memSingleBlock* queryPositionList;
    struct queryPosition* queryPositionEntry;
    struct codeword* oldCodeword;

    // Record word lookup structure size details
    wordLookupDFA_numCodes = numCodes;
    wordLookupDFA_wordLength = wordLength;

    // Start with no additional query positions
	wordLookupDFA_additionalQueryPositions = NULL;
    wordLookupDFA_numAdditionalQueryPositions = 0;

    // Maximum number of entry codewords and group codewords
    numWords = ceil(pow(numCodes, wordLength));
	numGroups = ceil(pow(numCodes, wordLength - 1));

	wordLookupDFA_numWords = numWords;
	wordLookupDFA_numGroups = numGroups;

    // Declare memory to construct initial words
	initialWords = (struct initialWord*)global_malloc(sizeof(struct initialWord)
                 * numWords);

    // Declare memory to hold list of neighbours for current query word
    if (encoding_alphabetType != encoding_nucleotide)
    	neighbours = (struct neighbour*)global_malloc(sizeof(struct neighbour) * numWords);

    // Stage 1 - Build initial word structure which contains, for each gram
    // a list of query positions that are high scoring

    // Iterate through every possible codeword
    codeword = 0;
    while (codeword < numWords)
    {
        // Initialize list of query positions as empty
		initialWords[codeword].queryPositions = NULL;
		initialWords[codeword].numQueryPositions = 0;
		initialWords[codeword].allocQueryPositions = 0;

        codeword++;
    }

    // Slide a word-sized window across the query
    queryPosition = 0;
    while (queryPosition < PSSMatrix.length - wordLength + 1)
    {
        // Get neighbours for the current query word
        numNeighbours = 0;
        wordLookupDFA_getNeighbours(PSSMatrix, queryPosition, &numNeighbours,
                                    neighbours);

        while (numNeighbours > 0)
        {
           	numNeighbours--;
            codeword = neighbours[numNeighbours].codeword;

            initialWord = initialWords + codeword;

            initialWord->numQueryPositions++;

            // Allocate memory to store additional query position if needed
            if (initialWord->numQueryPositions > initialWord->allocQueryPositions)
            {
                initialWord->allocQueryPositions
                    += constants_initialAllocCodewordQueryPositions;
                initialWord->queryPositions = (int2*)global_realloc(initialWord->queryPositions,
                    sizeof(int2) * initialWord->allocQueryPositions);
            }

            // Add query position to codeword (record position at END of word)
            initialWord->queryPositions[initialWord->numQueryPositions - 1] =
                queryPosition + wordLength - 1;

            // Record total number of query positions that will not fit into
            // a unsigned char and will need to be stored externally
            if (initialWord->numQueryPositions == 1)
            {
                numExternalPositions += 2;
            }
            else
            {
                numExternalPositions++;
            }
        }

        queryPosition++;
    }

    free(neighbours);

    // Contruct array of list of groups that share the same first query position
    groupsHash = (struct memSingleBlock**)global_malloc(sizeof(struct memSingleBlock*) * PSSMatrix.length);

    // Initialize each array cell to an empty list
    queryPosition = 0;
    while (queryPosition < PSSMatrix.length)
    {
		groupsHash[queryPosition] = memSingleBlock_initialize(sizeof(struct memSingleBlock*), 10);
    	queryPosition++;
	}

    // Go through each group
    groupCodeword = 0;
    while (groupCodeword < numGroups)
    {
    	// Find the first query position in the group
        code = 0;
        queryPosition = 0;
        while (code < numCodes && queryPosition == 0)
        {
            codeword = groupCodeword * numCodes + code;

            // Found an entry
            if (initialWords[codeword].numQueryPositions > 0)
            {
				queryPosition = initialWords[codeword].queryPositions[0];
            }

            code++;
        }

        // Add group's code to hash using that query position
        // (If no query positions, add to hash under the zero entry)
        groupList = memSingleBlock_initialize(sizeof(struct memSingleBlock*), 1);
        *((int4*)memSingleBlock_newEntry(groupList)) = groupCodeword;
		*((struct memSingleBlock**)memSingleBlock_newEntry(groupsHash[queryPosition])) = groupList;

		groupCodeword++;
    }

    // For each entry in the hash table
    queryPosition = 0;
    while (queryPosition < PSSMatrix.length)
    {
    	groupLists = (struct memSingleBlock**)groupsHash[queryPosition]->block;
    	// For each group entry A in that hash
        count = 0;
        while (count < groupsHash[queryPosition]->numEntries)
        {
			// For each group entry B after A in the hash
            count2 = count + 1;
			while(count2 < groupsHash[queryPosition]->numEntries)
            {
                // If both entries contain at least one group
                if (groupLists[count]->numEntries > 0 && groupLists[count2]->numEntries > 0)
                {
                    // Check if they have matching query positions
                    match = 1;
                    code = 0;
                    while (code < numCodes)
                    {
                        // For each word in the group
                        groupList1 = (int4*)(groupLists[count]->block);
                        groupList2 = (int4*)(groupLists[count2]->block);

                        codeword = groupList1[0] * numCodes + code;
                        codeword2 = groupList2[0] * numCodes + code;

                        // Check if the number of positions matches
                        if (initialWords[codeword].numQueryPositions !=
                            initialWords[codeword2].numQueryPositions)
                        {
                            // No match, break out
                            match = 0;
                            break;
                        }

                        // Move though the query positions and check for matches
                        queryPositionCount = 0;
                        while (queryPositionCount < initialWords[codeword].numQueryPositions)
                        {
                            // If they don't match at this position
                            if (initialWords[codeword].queryPositions[queryPositionCount] !=
                                initialWords[codeword2].queryPositions[queryPositionCount])
                            {
                                // No match, break out
                                match = 0;
                                break;
                            }
                            queryPositionCount++;
                        }

                        code++;
                    }

                    if (match)
                    {
                        // Move group in entry B to entry A
                        *((int4*)memSingleBlock_newEntry(groupLists[count])) = groupList2[0];

                        // Remove group from entry B
                        groupLists[count2]->numEntries = 0;
                    }
				}

            	count2++;
            }

        	count++;
        }

    	queryPosition++;
	}

    // Calculate frequencies of all possible groups of codes length one less than
    // the full word length and then arrange them so the most frequent have central positions
    aaFrequencyGroups = wordLookupDFA_calcFrequencyGroups(groupsHash, PSSMatrix,
                                                          numCodes, wordLength - 1);

    // Declare memory for addresses of group sized blocks
    wordLookupDFA_blockAddresses = (unsigned char**)global_malloc(sizeof(unsigned char*) * numGroups);

    // Declare memory for final words
	wordLookupDFA = (unsigned char*)global_malloc(sizeof(char) * numWords + sizeof(int2)
                  * numExternalPositions + sizeof(int2) * numExternalPositions);

    // Declare memory to store additional query positions
    // (we assume worst case scenario where all positions are addition positions)
    wordLookupDFA_additionalQueryPositions
    	= (uint2*)((char*)wordLookupDFA + sizeof(char) * numWords
          + sizeof(int2) * numExternalPositions);

	//wordLookupDFA_offsetFP = sizeof(char) * numWords + sizeof(int2) * numExternalPositions;
	wordLookupDFA_numExtPositions = numExternalPositions;

	currentBlock = wordLookupDFA;

    // Initialize query position list construction
    qPosList_initialize(numCodes);

    // Stage 2 - Construct the more compact, final wordLookupDFA structure that is arranged
    // with high frequency words in the middle (protein only)

    // For each collection of groups with the same set of query positions
    count = 0;
	while (count < wordLookupDFA_numBlocks)
    {
        groupList = aaFrequencyGroups[count].groups;
        groupList1 = (int4*)(groupList->block);

        // Use first group entry for constructing the block
        groupCodeword = groupList1[0];

        // First comes the external codes.
        startExternalPositions = (int2*)(currentBlock);
        endExternalPositions = startExternalPositions;

		qPosList_reset();

        // For each word in this block
        code = 0;
        while (code < numCodes)
        {
        	// Determine codeword
            codeword = groupCodeword * numCodes + code;
            queryPositionCount = initialWords[codeword].numQueryPositions;

            // Add list of query positions for this word to the qPosList structure
			qPosList_addList(initialWords[codeword].queryPositions,
                             initialWords[codeword].numQueryPositions, code);

			code++;
        }

        // Process the list of query positions to create an optimal set of lists
		qPosList_processLists();

        // For each query position list in order from shortest to longest
        while (qPosList_numQPosLists > 0)
        {
        	qPosList_numQPosLists--;
            queryPositionList = qPosList_qPosLists + qPosList_numQPosLists;

            // First check that there is room for the list in external area,
            // plus at least 3 int2s for each remaining code all within
            // the range of an unsigned char (256)
            if (endExternalPositions - startExternalPositions + queryPositionList->numEntries
                + 1 >= 256 - 3 * numCodes)
			{
                // If not, we need to store this word's query positions
                // at a completely external location

                // Copy values to additional positions address
                count2 = queryPositionList->numEntries;
                while (count2 > 0)
                {
                	count2--;
                	queryPositionEntry = memSingleBlock_getEntry(queryPositionList, count2);

					// Copy the query position
                    wordLookupDFA_additionalQueryPositions[
                    	wordLookupDFA_numAdditionalQueryPositions]
                        = queryPositionEntry->queryPosition;

                    // If a codeword list starts at this point
					while (queryPositionEntry->codewords != NULL)
                    {
                    	// Get each code
                    	code = queryPositionEntry->codewords->codeword;

                        // Add address information to external positions
                        // First flag to indicate additional positions
                        *endExternalPositions = 0;

                        // That's where the word->queryPosition will point to
                        externalPositions[code] = endExternalPositions;
                        endExternalPositions++;

                        // Then the offset into the array where they start
                        // (Record across two int2's to support very long queries)
                        *endExternalPositions = wordLookupDFA_numAdditionalQueryPositions / constants_max_int2;
                        endExternalPositions++;
                        *endExternalPositions = wordLookupDFA_numAdditionalQueryPositions % constants_max_int2;
                        endExternalPositions++;

                        // Free the current codeword before moving to the next
                        oldCodeword = queryPositionEntry->codewords;
                        queryPositionEntry->codewords = queryPositionEntry->codewords->next;
                        free(oldCodeword);
                    }

	                // Increment number of additional query positions
                    wordLookupDFA_numAdditionalQueryPositions++;
                }

                // Add terminating zero
                wordLookupDFA_additionalQueryPositions[
                    wordLookupDFA_numAdditionalQueryPositions] = 0;

                // Increment number of additional query positions
                wordLookupDFA_numAdditionalQueryPositions++;
            }
            else
            {
            	// Use external location for list of query positions

                // Copy values to external positions address
                count2 = queryPositionList->numEntries;
                while (count2 > 0)
                {
                	count2--;
                	queryPositionEntry = memSingleBlock_getEntry(queryPositionList, count2);

					// Copy the query position
                    *endExternalPositions = queryPositionEntry->queryPosition;

                    // If a codeword list starts at this point
					while (queryPositionEntry->codewords != NULL)
                    {
                    	// Get each code
                    	code = queryPositionEntry->codewords->codeword;

                        // The word will point to this external position
                        externalPositions[code] = endExternalPositions;

                        // Free the current codeword before moving to the next
                        oldCodeword = queryPositionEntry->codewords;
                        queryPositionEntry->codewords = queryPositionEntry->codewords->next;
                        free(oldCodeword);
					}

                    endExternalPositions++;
				}

                // Add zero terminal code
                *endExternalPositions = 0;
                endExternalPositions++;
            }
        }

        // Next comes the block structure
        currentPosition = (unsigned char*)endExternalPositions;

        // For each group, point to that block
        groupCount = 0;
        while (groupCount < groupList->numEntries)
        {
            wordLookupDFA_blockAddresses[groupList1[groupCount]] = currentPosition;
            groupCount++;
        }

        // Free the group list
        memSingleBlock_free(groupList);

        code = 0;
        while (code < numCodes)
        {
            codeword = groupCodeword * numCodes + code;

            queryPositionCount = initialWords[codeword].numQueryPositions;

            // No query positions, store NULL both 16-bit integers
            if (queryPositionCount == 0)
            {
                *currentPosition = 0;
            }
            // At least one query position
            else
            {
                // Give location of external query positions
                *currentPosition
                	= (unsigned char)(((char*)endExternalPositions -
                                       (char*)externalPositions[code]) / sizeof(int2));
            }

            // Free memory used by initial structure
            free(initialWords[codeword].queryPositions);

            code++;
            currentPosition++;
        }
        currentBlock = currentPosition;

        count++;
	}

    free(initialWords);
    free(aaFrequencyGroups);

    // Free the query position lists structure
    qPosList_free();

    // New stage 3 - construct groups lookup FSA structure

    // Declare memory for group states
    wordLookupDFA_groups = (struct group*)global_malloc(sizeof(struct group) * numGroups);

    // For each group
    groupCodeword = 0;
    while (groupCodeword < numGroups)
    {
        // Get block associated with current group
        currentBlock = wordLookupDFA_blockAddresses[groupCodeword];

		// Set pointer to block of words
        wordLookupDFA_groups[groupCodeword].nextWords = currentBlock;

        // Set pointer to next groups
		wordLookupDFA_groups[groupCodeword].nextGroups =
			wordLookupDFA_groups + ((groupCodeword * numCodes) % numGroups);

        groupCodeword++;
	}

	//Shucai
	//Transform DFP Lookup table to the format compatible with GPU
	wordLookupDFA_transformationForGPU(numCodes, wordLength);
}

//Shucai
//transform the current lookup table that one that can be used on GPU
void wordLookupDFA_transformationForGPU(int4 numCodes, int4 wordLength)
{
	int4 numWords, numGroups;
	int4 groupCodeword;

	//Maximum number of entry codewords and group codewords
	numWords = ceil(pow(numCodes, wordLength));
	numGroups = ceil(pow(numCodes, wordLength - 1));

	//Allocate wordLookupDFA_groupsFP buffer
	wordLookupDFA_groupsFP = (struct groupFP *)global_malloc(sizeof(struct groupFP) * numGroups);

	//Transform linklist table to an Array table
	groupCodeword = 0;
	while (groupCodeword < numGroups)
	{
		wordLookupDFA_groupsFP[groupCodeword].nextWords = (wordLookupDFA_groups[groupCodeword].nextWords - wordLookupDFA) / sizeof(char);
		wordLookupDFA_groupsFP[groupCodeword].nextGroups = (wordLookupDFA_groups[groupCodeword].nextGroups - wordLookupDFA_groups);
		groupCodeword++;
	}
	printf("Lookup table transformation completed...\n");
}

// Recursively find neighbours of the current codeword
void wordLookupDFA_findNeighbours(struct PSSMatrix PSSMatrix, int4 queryPosition, int4* numNeighbours,
                                  struct neighbour* neighbours)
{
	unsigned char* codes;
	unsigned char oldValue, tempValue;
    int4 position, oldScore, newScore, newCodeword, currentNeighbourCount;
	int4 currentScore, codeword;

	// For each neighbour in the growing list
    currentNeighbourCount = 0;
    while (currentNeighbourCount < *numNeighbours)
    {
    	codeword = neighbours[currentNeighbourCount].codeword;
        currentScore = neighbours[currentNeighbourCount].score;

        // Start at neighbour's current position
        position = neighbours[currentNeighbourCount].position;

        // Convert codeword to codes
        codes = wordLookupDFA_getCodes(codeword, wordLookupDFA_wordLength);

        // Move through each position not yet altered and try changing the code at that position
        while (position < wordLookupDFA_wordLength)
        {
            // Record old code at this position, and old score associated with it
            oldValue = codes[position];

            oldScore = PSSMatrix.matrix[queryPosition + position][oldValue];

            // For each code
            tempValue = 0;
            while (tempValue < wordLookupDFA_numCodes)
            {
                // If not the existing value at that position
                if (tempValue != oldValue)
                {
                    // Change the code
                    codes[position] = tempValue;

                    // Test if word is high scoring
                    newScore = PSSMatrix.matrix[queryPosition + position][tempValue];

                    if (currentScore - oldScore + newScore >= parameters_T)
                    {
                        // If so get the codeword
                        newCodeword = wordLookupDFA_getCodeword(codes, wordLookupDFA_wordLength);

                        neighbours[*numNeighbours].codeword = newCodeword;
                        neighbours[*numNeighbours].score = currentScore - oldScore + newScore;
                        neighbours[*numNeighbours].position = position + 1;
                        (*numNeighbours)++;
                    }
                }
                tempValue++;
            }

            // Change code back to old value before moving to next position
            codes[position] = oldValue;

            position++;
        }

        free(codes);

        currentNeighbourCount++;
	}
}

// Get all the neighbours for given query window
void wordLookupDFA_getNeighbours(struct PSSMatrix PSSMatrix, int4 queryPosition,
                                 int4* numNeighbours, struct neighbour* neighbours)
{
	int4 codeword, score = 0, count = queryPosition, containsWild = 0;

    // Get score for aligning the best match codes to the query window
	while (count < queryPosition + wordLookupDFA_wordLength)
    {
    	if (PSSMatrix.queryCodes[count] >= encoding_numRegularLetters)
			containsWild = 1;

        score += PSSMatrix.matrix[count][PSSMatrix.bestMatchCodes[count]];
        count++;
    }

	// If a word containing wildcards only consider nearest neighbour if high scoring
	if (!containsWild || score >= parameters_T)
    {
        // Convert query word codes to codeword
        codeword = wordLookupDFA_getCodeword(PSSMatrix.bestMatchCodes + queryPosition,
                                             wordLookupDFA_wordLength);

        // Automatically add the query word itself to list of neighbours
        neighbours[*numNeighbours].codeword = codeword;
        neighbours[*numNeighbours].score = score;
        neighbours[*numNeighbours].position = 0;
        (*numNeighbours)++;

        // Recursively find remaining neighbours
        wordLookupDFA_findNeighbours(PSSMatrix, queryPosition, numNeighbours, neighbours);
	}
}

int4 wordLookupDFA_compareFrequencyPair(const void* pair1, const void* pair2)
{
	const struct aaFrequencyGroup *p1, *p2;

	p1 = (struct aaFrequencyGroup*)pair1;
	p2 = (struct aaFrequencyGroup*)pair2;

	if (p1->frequency > p2->frequency)
	{
		return -1;
	}
	else if (p1->frequency < p2->frequency)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

// Given a total number of codes, and a word length, find all groups of codes that
// are wordLength in length and arrange them according to their background frequency
struct aaFrequencyGroup* wordLookupDFA_calcFrequencyGroups(struct memSingleBlock** groupsHash,
                              struct PSSMatrix PSSMatrix, int4 numCodes, int4 wordLength)
{
	struct aaFrequencyGroup* aaFrequencyGroups;
	struct aaFrequencyGroup* arrangedGroups;
    int4 count, middle, numGroups;
	unsigned char *codes;
    uint4 codeCount, blockCount = 0, queryPosition, groupCount;
	struct memSingleBlock **groupLists;
    int4* groupList1;
    float frequency;

    // Maximum number of group codewords
    numGroups = ceil(pow(numCodes, wordLength));

    // Construct amino acid frequencies array
    aaFrequencyGroups = (struct aaFrequencyGroup*)global_malloc(sizeof(struct aaFrequencyGroup)
                      * numGroups);
	arrangedGroups = (struct aaFrequencyGroup*)global_malloc(sizeof(struct aaFrequencyGroup)
                   * numGroups);

    // Iterate through each entry in the groups hash
    queryPosition = 0;
    while (queryPosition < PSSMatrix.length)
    {
    	groupLists = (struct memSingleBlock**)groupsHash[queryPosition]->block;

        // For each set of groups with identical query positions
        count = 0;
        while (count < groupsHash[queryPosition]->numEntries)
        {
        	// If there are groups at this entry
            if (groupLists[count]->numEntries > 0)
			{
            	wordLookupDFA_numBlocks++;
                groupList1 = (int4*)(groupLists[count]->block);

                aaFrequencyGroups[blockCount].frequency = 0;
                aaFrequencyGroups[blockCount].groups = groupLists[count];

                // For each group, determine the frequency and add to the total
                groupCount = 0;
                while (groupCount < groupLists[count]->numEntries)
                {
                	codes = wordLookupDFA_getCodes(groupList1[groupCount], wordLength);

                    // Determine the frequency for this group
	                frequency = 1;
                    codeCount = 0;
                    while (codeCount < wordLength)
                    {
                        frequency *= Robinson_prob[codes[codeCount]];

                        codeCount++;
                    }

                    free(codes);

                    // Add to the total
                    aaFrequencyGroups[blockCount].frequency += frequency;

                    groupCount++;
                }

                blockCount++;
            }
            else
            {
            	// Just free the entry, it's empty anyway
                memSingleBlock_free(groupLists[count]);
            }
        	count++;
        }

        // Free the current query position list
        memSingleBlock_free(groupsHash[queryPosition]);

    	queryPosition++;
	}

    // Free the list of query positions
	free(groupsHash);

    // Sort frequency pairs
	qsort(aaFrequencyGroups, blockCount, sizeof(struct aaFrequencyGroup),
	      wordLookupDFA_compareFrequencyPair);

    // Arrange into second array so that high frequency pairs have a central location
    // and lower frequency ones are on the outside
    middle = floor((float)(blockCount - 1) / 2);
    count = 0;
    while (count < blockCount)
    {
    	if (count % 2 == 0)
        {
        	arrangedGroups[middle - count / 2] = aaFrequencyGroups[count];
        }
        else
        {
        	arrangedGroups[middle + (count + 1) / 2] = aaFrequencyGroups[count];
        }

        count++;
    }

    free(aaFrequencyGroups);

	wordLookupDFA_numBlocks = blockCount;

    return arrangedGroups;
}

// Given a set of codes, construct and return the codeword
inline int4 wordLookupDFA_getCodeword(unsigned char* codes, int4 wordLength)
{
	int4 codeword;
    uint4 codeCount;

    codeword = 0;
    codeCount = 0;
    while (codeCount < wordLength)
    {
        codeword *= wordLookupDFA_numCodes;
        codeword += codes[codeCount];
        codeCount++;
    }

    return codeword;
}

// Given a code word, returns a list of codes
inline unsigned char* wordLookupDFA_getCodes(int4 codeword, int4 wordLength)
{
	unsigned char* codes;
	int4 count = wordLength;

    codes = (unsigned char*)global_malloc(sizeof(unsigned char) * wordLength);

	while (count > 0)
    {
    	count--;
    	codes[count] = codeword % wordLookupDFA_numCodes;
        codeword /= wordLookupDFA_numCodes;
    }

    return codes;
}

void wordLookupDFA_printCodes(int4 codeword, int4 wordLength)
{
	unsigned char* codes;
    int4 count = 0;

    codes = wordLookupDFA_getCodes(codeword, wordLength);

	while (count < wordLength)
    {
    	printf("%c", encoding_getLetter(codes[count]));
    	count++;
    }
    printf("\n");

    free(codes);
}

// Print the contents of the word lookup table
void wordLookupDFA_print()
{
    unsigned char *currentBlock;
	unsigned char codes[wordLookupDFA_wordLength + 1], code;
	unsigned char word;
	int4 groupCodeword;
	uint2* queryPositions;
    int4 totalQueryPositions = 0, totalEmptySlots = 0;
	int4 count;
    int4 numWords, numGroups;
    int4 numCodes, wordLength;
    uint4 codeCount;

	numCodes = wordLookupDFA_numCodes;
    wordLength = wordLookupDFA_wordLength;

    // Determine number of codewords and number of groups
    numGroups = ceil(pow(numCodes, wordLength - 1));
    numWords = ceil(pow(numCodes, wordLength));

    // Initialize word to first position
    codeCount = 0;
    while (codeCount <= wordLength - 1)
    {
    	codes[codeCount] = 0;
        codeCount++;
    }

    // Iterate - For each possible group
	while (!codes[wordLength - 1])
    {
    	// Construct the codeword for array of codes
		groupCodeword = wordLookupDFA_getCodeword(codes, wordLength - 1);

        // Get current block
        currentBlock = wordLookupDFA_blockAddresses[groupCodeword];

        // For each word in the group
        code = 0;
        while (code < numCodes)
        {
            // Get word in the block
            word = currentBlock[code];

            if (encoding_alphabetType == encoding_protein || word != 0)
            {
                //  Print the word
                codeCount = 0;
                while (codeCount < wordLength - 1)
                {
                    printf("%c", encoding_getLetter(codes[codeCount]));
                    codeCount++;
                }
                printf("%c", encoding_getLetter(code));

                printf(" QueryPositions ");
                fflush(stdout);

                printf("[%d]: ", word);

                if (word == 0)
                {
                    // No query positions
                    totalEmptySlots++;
                }
                else
                {
                    // At least one position at an external address
                    queryPositions = ((uint2*)currentBlock) - word;

                    // If the zero flag is stored at the first query position
                    if (!*queryPositions)
                    {
                    	// Go to an outside address for additional positions
						queryPositions = wordLookupDFA_additionalQueryPositions
                                       + *(queryPositions + 1);
                    }

                    count = 0;
                    while (queryPositions[count] != 0)
                    {
                        printf("%d ", queryPositions[count] - 1);
                        totalQueryPositions++;
                        fflush(stdout);
                        count++;
                    }
                }

                printf("\n");
            }

            code++;
		}

        // Move to next word
		codes[0]++;
    	codeCount = 0;
        while (codeCount < wordLength - 1)
        {
			if (codes[codeCount] >= numCodes)
            {
            	codes[codeCount] = 0;
                codes[codeCount + 1]++;
            }
            else
            {
            	break;
            }
            codeCount++;
        }
    }

    printf("Total query positions=%d\n", totalQueryPositions);
    printf("Empty slots=%d/%d (%d%%)\n", totalEmptySlots, numWords,
                                         totalEmptySlots * 100 / numWords);
	printf("Number of blocks/groups=%d/%d\n", wordLookupDFA_numBlocks, numGroups);
}

// Free memory used by the word lookup table
void wordLookupDFA_free()
{
	free(wordLookupDFA_groups);
    free(wordLookupDFA);
    free(wordLookupDFA_blockAddresses);
	//Shucai
	free(wordLookupDFA_groupsFP);
}
