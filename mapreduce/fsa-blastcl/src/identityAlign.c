// identityAlign.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform Smith-Waterman alignment using 0-value speedup used by
// Phil Green's implementation. Returns score and end-points only.

#include "blast.h"

int4* identityAlign_bestRow = NULL;
int4 identityAlign_rowSizes = 0;

#define maximum(a,b) ((a > b) ? a : b)

void identityAlign_print(int4 bandStart, int4 bandEnd, int4 seqEnd);

// Perform banded smith-waterman dynamic programming on the given query and subject sequences
// where +1 is rewarded for a match, no penalty for mismatch or gap.
// seq1 must be shorter than seq2 and alignment will stop if score cannot reach
// length1 - maxMismatch
int4 identityAlign_score(unsigned char* seq1, int4 length1,
                         unsigned char* seq2, int4 length2, int4 relativeOffset, uint4 targetScore)
{
	unsigned char *position1, *bestPosition1, *end1;
	unsigned char *position2, *bestPosition2, *end2, *bandStart, *bandEnd;
    int4* bestRow;
	int4 oldBest, match, previousOldBest;
	int4 previousBest;
	int4 bestScore = 0, rowsRemaining = length1;

//    printf("lengths=%d,%d\n", length1, length2);

    bandStart = seq2 - relativeOffset - 10;
    bandEnd = seq2 - relativeOffset + 10;

    end1 = seq1 + length1;
    end2 = seq2 + length2;

	// Declare processing rows for storing match, insert-subject and insert-query values
	// If current malloced rows aren't big enough
	if (length2 >= identityAlign_rowSizes)
	{
		// Free existing rows
		free(identityAlign_bestRow);
		// Set size to double current needed length
		identityAlign_rowSizes = length2 + 1;
		// Malloc new rows
		identityAlign_bestRow = (int4*)global_malloc(sizeof(int4) * identityAlign_rowSizes);
	}

	position1 = seq1 - 1;
	position2 = seq2 - 1;

	// Initialize rows
	bestRow = identityAlign_bestRow;

	// -----FIRST ROW-----
	// For each cell in the top row, scanning from left-to-right
	while (position2 < end2)
	{
		// All values are zero
		*bestRow = 0;

		bestRow++; position2++;
	}

	position1++;

	// -----REMAINING ROWS-----
	while (position1 < end1)
	{
        bandStart++;
        bandEnd++;

		position2 = maximum(seq2 - 1, bandStart);
		bestRow = identityAlign_bestRow + (position2 - seq2) + 1;

        if (position2 < end2 && position2 < bandEnd)
        {
            // -----FAR LEFT CELL-----
            // Set zero values
            *bestRow = 0;
            previousOldBest = 0;
            previousBest = 0;
		}

		position2++; bestRow++;

		// -----REMAINING CELLS-----
		while (position2 < end2 && position2 < bandEnd)
		{
			// Remember old Best value (for cell below this one)
			oldBest = *bestRow;

			// Calculate new M value
            if (*position1 == *position2)
            	match = previousOldBest + 1;
            else
            	match = previousOldBest;

			previousOldBest = oldBest;

            // Determine the best of M and Ix
            *bestRow = maximum(maximum(match, oldBest), previousBest);

            previousBest = *bestRow;

            // If this is the best-yet scoring cell
            if (match > bestScore)
            {
                // Update best start cell data
                bestScore = match;
            }

            position2++; bestRow++;
		}

        position1++;
		rowsRemaining--;

        if (rowsRemaining + bestScore < targetScore)
        {
//        	printf("Stopped remaining %d/%d\n", rowsRemaining, length1);
        	return 0;
		}

//        identityAlign_print(bandStart - seq2, bandEnd - seq2, end2 - seq2);
//		if (dloc==72249)
//			smithWatermanTraceback_print(identityAlign_bestRow, subjectLength);
	}

//    printf("Finished %d\n", bestScore);
	return bestScore;
}

void identityAlign_print(int4 bandStart, int4 bandEnd, int4 seqEnd)
{
    int4 count = 0;
//	printf("[%d,%d,%d]", bandStart, bandEnd, seqEnd); fflush(stdout);

    while (count < bandStart)
    {
		printf("    ");
    	count++;
    }

    while (count < bandEnd && count < seqEnd)
    {
		printf("%3d ", identityAlign_bestRow[count]);
        count++;
    }

    printf("\n");
}

void identityAlign_free()
{
	free(identityAlign_bestRow);
}
