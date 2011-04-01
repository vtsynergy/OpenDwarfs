// nucleotideLookup.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for creating a nucleotide word lookup table

#include "blast.h"

int2* nucleotideLookup = NULL;
uint2* nucleotideLookup_additionalPositions = NULL;
uint4 nucleotideLookup_numAdditionalPositions;
uint4 *nucleotideLookup_bitLookup = NULL, nucleotideLookup_mask;
int4 nucleotideLookup_packedWordLength;

int4* nucleotideLookup_large = NULL;
uint4* nucleotideLookup_additionalPositions_large = NULL;
char nucleotideLookup_largeTable = 0;

struct initialWord
{
	int4 numQueryPositions;
	int4 allocQueryPositions;
    int4* queryPositions;
};

// Build the nucleotide lookup table
void nucleotideLookup_build(struct PSSMatrix PSSMatrix, int4 packedWordLength)
{
	uint4 codeword, numEntries, wordLength, byteNumber, queryPosition;
    uint4 numBitWords, numAdditionalPositions = 1;
	struct initialWord *initialLookup, *initialWord;
    int4 count;

    nucleotideLookup_numAdditionalPositions = 1;

    // Number of entries in the table
    numEntries = ceil(pow(256, packedWordLength));
    wordLength = packedWordLength * 4;
	nucleotideLookup_packedWordLength = packedWordLength;

    // Declare memory for initial lookup table
    initialLookup = (struct initialWord*)global_malloc(sizeof(struct initialWord) * numEntries);

    // Iterate through every possible codeword
    codeword = 0;
    while (codeword < numEntries)
    {
        // Initialize list of query positions as empty
		initialLookup[codeword].queryPositions = NULL;
		initialLookup[codeword].numQueryPositions = 0;
		initialLookup[codeword].allocQueryPositions = 0;

        codeword++;
    }

    // Slide a word-sized window across the query
    queryPosition = 0;
    while (queryPosition < PSSMatrix.length - wordLength + 1)
    {
    	codeword = 0;

        // Don't include words that cross the strand boundry
        if (!(queryPosition < PSSMatrix.strandLength &&
              queryPosition >= PSSMatrix.strandLength - wordLength + 1))
		{
            // For each set of four letters in the word
            byteNumber = 0;
            while (byteNumber < packedWordLength)
            {
                // Packed them and add to the codeword
                codeword = codeword << 8;
                codeword |= (int4)encoding_bytePack(PSSMatrix.queryCodes + queryPosition + byteNumber * 4);
                byteNumber++;
            }

            initialWord = initialLookup + codeword;
            initialWord->numQueryPositions++;

            // Allocate memory to store additional query position if needed
            if (initialWord->numQueryPositions > initialWord->allocQueryPositions)
            {
                initialWord->allocQueryPositions
                    += constants_initialAllocCodewordQueryPositions;
                initialWord->queryPositions = (uint4*)global_realloc(initialWord->queryPositions,
                    sizeof(uint4) * initialWord->allocQueryPositions);
            }

            // Add query position to codeword (record position at END of word)
            initialWord->queryPositions[initialWord->numQueryPositions - 1] = queryPosition + wordLength;

            if (initialWord->numQueryPositions == 2)
            {
            	numAdditionalPositions += 3;
            }
            else if (initialWord->numQueryPositions > 2)
            {
            	numAdditionalPositions++;
            }
        }

        queryPosition++;
    }

    // For long queries we need to use the larger 32-bit lookup table
    if (numAdditionalPositions > constants_max_int2 || PSSMatrix.length > constants_max_int2)
    {
    	nucleotideLookup_largeTable = 1;

        // Declare memory for table
        nucleotideLookup_large = (int4*)global_malloc(sizeof(int4) * numEntries);

        // Declare memory for additional query positions
        nucleotideLookup_additionalPositions_large = (uint4*)global_malloc(sizeof(uint4)
                                                   * numAdditionalPositions);
    }
    // Shorter queries we use the smaller 16-bit lookup table
    else
    {
		nucleotideLookup_largeTable = 0;

        // Declare memory for table
        nucleotideLookup = (int2*)global_malloc(sizeof(int2) * numEntries);

        // Declare memory for additional query positions
        nucleotideLookup_additionalPositions = (uint2*)global_malloc(sizeof(uint2)
                                             * numAdditionalPositions);
	}

	// Declare memory for bit-lookup table
    numBitWords = numEntries / 32;
	nucleotideLookup_bitLookup = (uint4*)global_malloc(sizeof(int4) * numBitWords);

    // Create a mask for the bit-lookup table
    nucleotideLookup_mask = 31;

    // For each word
    codeword = 0;
    while (codeword < numEntries)
    {
    	initialWord = initialLookup + codeword;

        // Set bit-lookup bit
		if (initialWord->numQueryPositions > 0)
        {
			nucleotideLookup_bitLookup[codeword >> 5] |= 1 << (codeword & nucleotideLookup_mask);
        }

        // No query positions
    	if (initialWord->numQueryPositions == 0)
        {
        	if (nucleotideLookup_largeTable)
				nucleotideLookup_large[codeword] = 0;
            else
				nucleotideLookup[codeword] = 0;
        }
        // One query position
    	else if (initialWord->numQueryPositions == 1)
        {
        	// Store in table
        	if (nucleotideLookup_largeTable)
				nucleotideLookup_large[codeword] = initialWord->queryPositions[0];
			else
				nucleotideLookup[codeword] = initialWord->queryPositions[0];

            free(initialWord->queryPositions);
        }
        // Multiple query positions
		else
        {
        	// Mark start of query positions
        	if (nucleotideLookup_largeTable)
            {
                nucleotideLookup_large[codeword] = -nucleotideLookup_numAdditionalPositions;

                // Copy to additional query position slots
                count = 0;
                while (count < initialWord->numQueryPositions)
                {
                    nucleotideLookup_additionalPositions_large[nucleotideLookup_numAdditionalPositions]
                        = initialWord->queryPositions[count];
                    count++;
                    nucleotideLookup_numAdditionalPositions++;
                }

                // Add null terminator
                nucleotideLookup_additionalPositions_large[nucleotideLookup_numAdditionalPositions] = 0;
                nucleotideLookup_numAdditionalPositions++;
			}
			else
            {
                nucleotideLookup[codeword] = -nucleotideLookup_numAdditionalPositions;

                // Copy to additional query position slots
                count = 0;
                while (count < initialWord->numQueryPositions)
                {
                    nucleotideLookup_additionalPositions[nucleotideLookup_numAdditionalPositions]
                        = initialWord->queryPositions[count];
                    count++;
                    nucleotideLookup_numAdditionalPositions++;
                }

                // Add null terminator
                nucleotideLookup_additionalPositions[nucleotideLookup_numAdditionalPositions] = 0;
                nucleotideLookup_numAdditionalPositions++;
			}

            free(initialWord->queryPositions);
        }

        codeword++;
    }

    free(initialLookup);
}

// Print the lookup table
void nucleotideLookup_print()
{
	int4 count, codeword, numEntries;

    numEntries = ceil(pow(256, nucleotideLookup_packedWordLength));

    // Print the table
    codeword = 0;
    while (codeword < numEntries)
    {
    	// Print the codeword
        if (nucleotideLookup[codeword])
        {
            count = nucleotideLookup_packedWordLength * 4;
            while (count > 0)
            {
                count--;
            	printf("%c", encoding_getLetter((codeword >> (count * 2)) & 0x3));
            }
		}

        if (nucleotideLookup[codeword] > 0)
        {
	        printf(": %d\n", nucleotideLookup[codeword]);
        }
        else if (nucleotideLookup[codeword] < 0)
        {
        	count = -nucleotideLookup[codeword];

	        printf(": ");

            while (nucleotideLookup_additionalPositions[count] != 0)
            {
            	printf("%d ", nucleotideLookup_additionalPositions[count]);
            	count++;
            }

            printf("\n");
        }
        codeword++;
    }
}

// Free the table
void nucleotideLookup_free()
{
    free(nucleotideLookup);
    free(nucleotideLookup_bitLookup);
    free(nucleotideLookup_additionalPositions);
    nucleotideLookup = NULL;
    nucleotideLookup_bitLookup = NULL;
    nucleotideLookup_additionalPositions = NULL;
}

