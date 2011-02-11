// smithWatermanScoring.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform Smith-Waterman alignment using 0-value speedup used by
// Phil Green's implementation. Returns score and end-points only.

#include "blast.h"

int4* smithWatermanScoring_bestRow = NULL;
int4* smithWatermanScoring_insertQrow = NULL;
int4 smithWatermanScoring_rowSizes = 0;

// Prototypes
struct dpResults smithWatermanScoring_dynamicProgramming(unsigned char* subject,
                                  struct PSSMatrix PSSMatrix, int4 subjectLength);
#define maximum(a,b) ((a > b) ? a : b)

// Perform smith-waterman dynamic programming on the given query and subject sequences
struct dpResults smithWatermanScoring_score(struct PSSMatrix PSSMatrix, int4 subjectSize,
                                            unsigned char* subject)
{
	struct dpResults dpResults, dpResults2;
	struct PSSMatrix choppedPSSMatrix;

    // If there are two strands
    if (PSSMatrix.strandLength != PSSMatrix.length)
    {
    	// Perform dynamic programming on first strand
        choppedPSSMatrix = PSSMatrix;
        choppedPSSMatrix.length = PSSMatrix.strandLength;
		dpResults = smithWatermanScoring_dynamicProgramming(subject, choppedPSSMatrix, subjectSize);

        // Perform dynamic programming on second strand
        choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, PSSMatrix.strandLength);
		dpResults2 = smithWatermanScoring_dynamicProgramming(subject, choppedPSSMatrix, subjectSize);
		dpResults2.best.queryOffset += PSSMatrix.strandLength;

        if (dpResults.bestScore > dpResults2.bestScore)
        	return dpResults;
        else
        {
//        	printf("[%d:%d,%d]", dpResults2.bestScore, dpResults2.best.queryOffset,
//                                 dpResults2.best.subjectOffset);
        	return dpResults2;
		}
    }
	else
    {
        // Perform dynamic programming
        dpResults = smithWatermanScoring_dynamicProgramming(subject, PSSMatrix, subjectSize);

//    	printf("Best=%d,%d Score=%d dloc=%d\n", dpResults.best.queryOffset, dpResults.best.subjectOffset,
//    	                                        dpResults.bestScore, blast_dloc); fflush(stdout);

        return dpResults;
	}
}

// Reverse the query and subject sequences, and perform dynamic programming to find the
// start of the optimal alignment, rather than the end
struct dpResults smithWatermanScoring_scoreReverse(struct PSSMatrix PSSMatrix, int4 subjectSize,
                                                   unsigned char* subject, struct coordinate end)
{
	unsigned char* reversedSubject;
	uint4 subjectPosition, strandOffset = 0;
    struct PSSMatrix reversedPSSMatrix;
	struct dpResults dpResults;

	end.queryOffset++;
	end.subjectOffset++;
	
    // If alignment is in the second strand
    if (end.queryOffset > PSSMatrix.strandLength)
    {
    	// Remove first strand
        strandOffset = PSSMatrix.strandLength;
        end.queryOffset -= PSSMatrix.strandLength;
        PSSMatrix = PSSMatrix_chop(PSSMatrix, PSSMatrix.strandLength);
	}

    // Remove the end of the PSSM not containing the alignment
    PSSMatrix.strandLength = PSSMatrix.length = end.queryOffset;

    // Shorten the subject to region containing the alignment
    subjectSize = end.subjectOffset;

    // Make a reversed copy of the subject region containing the alignment
    reversedSubject = (unsigned char*)global_malloc(sizeof(char) * subjectSize);
	subjectPosition = 0;
    while (subjectPosition < subjectSize)
    {
    	reversedSubject[subjectPosition] = subject[subjectSize - subjectPosition - 1];
    	subjectPosition++;
    }

    // Reverse the PSSM
    reversedPSSMatrix = PSSMatrix_reverse(PSSMatrix);

    // Perform smith waterman on reversed sequences
	dpResults = smithWatermanScoring_score(reversedPSSMatrix, subjectSize, reversedSubject);

    dpResults.best.queryOffset = PSSMatrix.length - dpResults.best.queryOffset - 1;
    dpResults.best.subjectOffset = subjectSize - dpResults.best.subjectOffset - 1;
	dpResults.best.queryOffset += strandOffset;

    PSSMatrix_freeCopy(reversedPSSMatrix);
    free(reversedSubject);

    return dpResults;
}

// Perform dynamic programming
struct dpResults smithWatermanScoring_dynamicProgramming(unsigned char* subject,
                                  struct PSSMatrix PSSMatrix, int4 subjectLength)
{
	int2 **queryPosition, **bestQueryPosition, **queryEnd;
	int2* matrixColumn;
	unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
	int4 bestScore = 0;
	int4 *bestRow, *insertQrow, insertS;
	int4 oldBest, match, previousOldBest;
	int4 previousBest;
	int4 queryLength;
	struct dpResults dpResults;
	int4 minimumForOpenGap;

	queryLength = PSSMatrix.length;
	subjectEnd = subject + subjectLength;
	queryEnd = PSSMatrix.matrix + queryLength;

	minimumForOpenGap = parameters_openGap;

	// Declare processing rows for storing match, insert-subject and insert-query values
	// If current malloced rows aren't big enough
	if (subjectLength >= smithWatermanScoring_rowSizes)
	{
		// Free existing rows
		free(smithWatermanScoring_bestRow);
		free(smithWatermanScoring_insertQrow);
		// Set size to double current needed length
		smithWatermanScoring_rowSizes = subjectLength + 1;
		// Malloc new rows
		smithWatermanScoring_bestRow = (int4*)global_malloc(sizeof(int4) * smithWatermanScoring_rowSizes);
		smithWatermanScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * smithWatermanScoring_rowSizes);
	}

	bestSubjectPosition = subject;
	subjectPosition = subject - 1;
	bestQueryPosition = PSSMatrix.matrix;
	queryPosition = PSSMatrix.matrix - 1;

	// Initialize rows
	bestRow = smithWatermanScoring_bestRow;
	insertQrow = smithWatermanScoring_insertQrow;

	// -----FIRST ROW-----

	// For each cell in the top row, scanning from left-to-right
	while (subjectPosition < subjectEnd)
	{
		// All values are zero
		*bestRow = *insertQrow = 0;

		bestRow++; insertQrow++; subjectPosition++;
	}

	queryPosition++;

	// -----REMAINING ROWS-----
	while (queryPosition < queryEnd)
	{
		subjectPosition = subject - 1;

		// Using next column of query matrix
		matrixColumn = *queryPosition;

		bestRow = smithWatermanScoring_bestRow;
		insertQrow = smithWatermanScoring_insertQrow;

		// -----FAR LEFT CELL-----

		// Set zero values
		previousOldBest = 0;
		*bestRow = *insertQrow = 0;
		match = insertS = previousBest = 0;

		subjectPosition++; bestRow++; insertQrow++;

		// -----REMAINING CELLS-----
start:
		// Loop 1: where insertS is zero
		while (subjectPosition < subjectEnd)
		{
			// Remember old Best value (for cell below this one)
			oldBest = *bestRow;

			// Calculate new M value
			match = maximum(matrixColumn[*subjectPosition] + previousOldBest, 0);
			previousOldBest = oldBest;

            // Match is not large enough to trigger opening a gap
            if (match <= minimumForOpenGap)
            {
                if (*insertQrow <= 0)
                {
                    *bestRow = match;
                }
                else
                {
                    *bestRow = maximum(*insertQrow, match);
                    *insertQrow -= parameters_extendGap;
                }
            }
            // Match is large enough to open a gap
            else
            {
                if (*insertQrow <= 0)
                {
                    // M is the only score computed, hence the best
                    *bestRow = match;

                    // Calculate new Iy
                    *insertQrow = insertS = match - parameters_openGap;
                }
                else
                {
                    // Determine the best of M and Ix
                    *bestRow = maximum(match, *insertQrow);

                    // Calculate new Ix
                    *insertQrow = maximum(match - parameters_openGap,
                                          *insertQrow - parameters_extendGap);

                    // Calculate new Iy
                    insertS = match - parameters_openGap;
                }

                // insertS > 0 so go to the second loop
                break;
            }

            // If this is the best-yet scoring cell
            if (match >= bestScore)
            {
                // Update best start cell data
                bestScore = match;
                bestQueryPosition = queryPosition;
                bestSubjectPosition = subjectPosition;
            }

            subjectPosition++; bestRow++; insertQrow++;
		}

        // Perform best-score-yet check before we move on
        if (match >= bestScore)
        {
            bestScore = match;
            bestQueryPosition = queryPosition;
            bestSubjectPosition = subjectPosition;
        }

        subjectPosition++; bestRow++; insertQrow++;

		// Loop 2: where insertS > 0
		while (subjectPosition < subjectEnd)
		{
			// Remember old Best value (for cell below this one)
			oldBest = *bestRow;

			// Calculate new M value
			match = maximum(matrixColumn[*subjectPosition] + previousOldBest, 0);
			previousOldBest = oldBest;

            // Match is not large enough to trigger opening a gap
            if (match <= minimumForOpenGap)
            {
                if (*insertQrow <= 0)
                {
                    *bestRow = maximum(insertS, match);
                }
                else
                {
                    *bestRow = maximum(match, maximum(insertS, *insertQrow));
                    *insertQrow -= parameters_extendGap;
                }

                // If S will not be above zero
            	if (insertS <= parameters_extendGap)
                {
                    // Go back to loop 1
                    break;
                }
				else
                {
                    insertS -= parameters_extendGap;
				}
            }
            // Match is large enough to open a gap
            else
            {
                if (*insertQrow <= 0)
                {
                    // Determine the best of M and Iy
                    *bestRow = maximum(match, insertS);

                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);

                    // Calculate new Ix
                    *insertQrow = match - parameters_openGap;
                }
                else
                {
                    // Determine the best of M, Ix and Iy
                    *bestRow = maximum(maximum(match, insertS), *insertQrow);

                    // Set new Ix value
                    *insertQrow = maximum(match - parameters_openGap,
                                          *insertQrow - parameters_extendGap);
                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);
                }
            }
			// If this is the best-yet scoring cell
			if (match >= bestScore)
			{
				// Update best start cell data
				bestScore = match;
                bestQueryPosition = queryPosition;
                bestSubjectPosition = subjectPosition;
			}

			subjectPosition++; bestRow++; insertQrow++;
		}

		subjectPosition++; bestRow++; insertQrow++;

        if (subjectPosition < subjectEnd)
        {
            goto start;
        }

        queryPosition++;

//		if (dloc==72249)
//			smithWatermanTraceback_print(smithWatermanScoring_bestRow, subjectLength);
	}

	dpResults.bestScore = bestScore;
    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.traceback = NULL;

//    	printf("Best=%d,%d Score=%d dloc=%d\n", dpResults.best.queryOffset, dpResults.best.subjectOffset,
//    	                                        dpResults.bestScore, blast_dloc); fflush(stdout);

	return dpResults;
}

