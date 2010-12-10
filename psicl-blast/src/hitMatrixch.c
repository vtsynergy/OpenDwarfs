// hitMatrix.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Array used by BLAST to detect two hits on the same diagonal, or in one-hit mode
// to detect if a hit has been covered by an existing ungapped extension

#include "blast.h"

// Initialize the hitMatrix by declaring memory for the maximum number
// of diagonals required by a subject sequence
void hitMatrix_initialize(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress)
{
	int4 numDiagonals;
	unsigned char **minOffset, **offset, **maxOffset;
	//Shucai
	uint4 *minOffsetFP, *offsetFP, *maxOffsetFP;

    // Use more memory efficient but slower hit matrix for nucleotide
    if (encoding_alphabetType == encoding_nucleotide)
    {
        // Calculate number of diagonals that will be required during search
        numDiagonals = 1;
        while (numDiagonals < queryLength + parameters_wordSize)
        {
            numDiagonals <<= 1;
        }

        // Construct mask
        hitMatrix_diagonalMask = numDiagonals - 1;

        // Declare memory for diagonal slots
        hitMatrix_furthest = (unsigned char**)global_malloc(sizeof(unsigned char*) * numDiagonals);
        minOffset = hitMatrix_furthest;
	}
    // Use less memory efficient but faster hit matrix for protein
    else
    {
        // Maximum number of diagonals that will be required during search
        numDiagonals = queryLength + maximumSubjectLength - parameters_wordSize + 1;
        minOffset = (unsigned char**)global_malloc(sizeof(unsigned char*) * numDiagonals);
		
		//Shucai
		minOffsetFP = (uint4 *)global_malloc(sizeof(uint4) * numDiagonals);

        // Advance array pointer to allow offset values ranging from
        // -queryLength to subjectLength - wordSize
        hitMatrix_furthest = minOffset + queryLength;
		//Shucai
		hitMatrix_furthestFP = minOffsetFP + queryLength;
    }

	// Record query length
	hitMatrix_queryLength = queryLength;

	// Start from smallest possible offset value and iterate through to largest
	offset = minOffset;
	maxOffset = minOffset + numDiagonals;

	// For each diagonal, reset furthest to address at start of file
	while (offset < maxOffset)
	{
		*offset = startAddress;
		offset++;
	}

	//Shucai
	offsetFP = minOffsetFP;
	maxOffsetFP = minOffsetFP + numDiagonals;
	while (offsetFP < maxOffsetFP)
	{
		*offsetFP = 0;
		offsetFP++;
	}
}

// Reinitialize the hit matrix so all values point to start of new file
void hitMatrix_reinitialize(int4 queryLength, int4 maximumSubjectLength, unsigned char* startAddress)
{
	int4 numDiagonals;
	unsigned char **minOffset, **offset, **maxOffset;
	//Shucai
	uint4 *minOffsetFP, *offsetFP, *maxOffsetFP;

    if (encoding_alphabetType == encoding_nucleotide)
    {
        // Calculate number of diagonals that will be required during search
        numDiagonals = 1;
        while (numDiagonals < queryLength + parameters_wordSize)
        {
            numDiagonals <<= 1;
        }

        minOffset = hitMatrix_furthest;
	}
    // Use less memory efficient but faster hit matrix for protein
    else
    {
        // Maximum number of diagonals that will be required during search
        numDiagonals = queryLength + maximumSubjectLength - parameters_wordSize + 1;

        minOffset = hitMatrix_furthest - queryLength;

		//Shucai
		minOffsetFP = hitMatrix_furthestFP - queryLength;
    }

	// Start from smallest possible offset value and iterate through to largest
	offset = minOffset;
	maxOffset = minOffset + numDiagonals;

	// For each diagonal, reset furthest to address at start of file
	while (offset < maxOffset)
	{
		*offset = startAddress;
		offset++;
	}

	//Shucai
	offsetFP = minOffsetFP;
	maxOffsetFP = minOffsetFP + numDiagonals;
	while (offsetFP < maxOffsetFP)
	{
		*offsetFP = 0;
		offsetFP++;
	}
}

// Free the matrix
void hitMatrix_free()
{
    if (encoding_alphabetType == encoding_protein)
		hitMatrix_furthest -= hitMatrix_queryLength;
	free(hitMatrix_furthest);

	//Shucai
	if (encoding_alphabetType == encoding_protein)
	{
		hitMatrix_furthestFP -= hitMatrix_queryLength;
	}
	free(hitMatrix_furthestFP);
}
