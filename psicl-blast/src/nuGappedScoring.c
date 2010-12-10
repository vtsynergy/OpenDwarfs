// gappedScoring.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform traditional gapped alignment WITH restricted insertion
// recording the score and alignment end-points only

#include "blast.h"

int4* nuGappedScoring_bestRow = NULL;
int4* nuGappedScoring_insertQrow = NULL;
int4 nuGappedScoring_rowSizes = 0;

#define maximum(a,b) ((a > b) ? a : b)

// Prototypes
void nuGappedScoring_printBeforeRow(int4* row, int2** query, int2** rowDropoff, int2** columnDropoff);
void nuGappedScoring_printAfterRow(int4* row, int2** query, int2** rowDropoff, int2** columnDropoff);
struct dpResults nuGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                              struct coordinate seed, int4 dropoff);
struct dpResults nuGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                             int4 dropoff, int4 subjectLength);

// Perform gapped alinment with restricted insertion
int4 nuGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                           int4 subjectSize, unsigned char* subject, int4 dropoff)
{
    struct coordinate seed;
    unsigned char *choppedSubject;
    struct PSSMatrix choppedPSSMatrix;
    int4 choppedSubjectSize;
    struct dpResults beforeDpResults, afterDpResults;
    int4 strandOffset = 0;

    seed = ungappedExtension->seed;
//    printf("Seed=%d,%d dloc=%d\n", seed.queryOffset, seed.subjectOffset, blast_dloc);
//	ungappedExtension_print(ungappedExtension); fflush(stdout);

    if (seed.queryOffset > PSSMatrix.strandLength)
    {
	    // If query position is in the second strand, remove first strand from PSSM
        strandOffset = PSSMatrix.strandLength;
		seed.queryOffset -= PSSMatrix.strandLength;
		PSSMatrix = PSSMatrix_chop(PSSMatrix, PSSMatrix.strandLength);
    }
    else
    {
    	// Otherwise remove second strand
    	PSSMatrix.length = PSSMatrix.strandLength;
    }

    // Adjust to point to start of byte
    seed.queryOffset -= 2;
    seed.subjectOffset -= 2;

    // Perform dynamic programming for points before the seed
    beforeDpResults = nuGappedScoring_dpBeforeSeed(subject, PSSMatrix, seed, dropoff);

    // Chop the start off the query and subject so they begin at the seed
    choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
    choppedSubject = subject + (seed.subjectOffset / 4);
    choppedSubjectSize = subjectSize - seed.subjectOffset;

    // Perform dynamic programming for points after the seed
    afterDpResults = nuGappedScoring_dpAfterSeed(choppedSubject, choppedPSSMatrix,
                                                 dropoff, choppedSubjectSize);

	// Convert best endpoints from bytepacked to unpacked offsets
	afterDpResults.best.subjectOffset *= 4;
	beforeDpResults.best.subjectOffset *= 4;

    // Re-adjust result change due to chopping subject/query and strand adjustment
    afterDpResults.best.queryOffset += seed.queryOffset + strandOffset;
    afterDpResults.best.subjectOffset += seed.subjectOffset;
    beforeDpResults.best.queryOffset += strandOffset;

    // Associate best scoring start and end points with the ungapped extension
    ungappedExtension->start = beforeDpResults.best;
    ungappedExtension->end = afterDpResults.best;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        printf("norm[%d,%d,%d](seed=%d,%d)\n", beforeDpResults.bestScore, afterDpResults.bestScore,
        	choppedPSSMatrix.matrix[0][encoding_extractBase(choppedSubject[0],seed.subjectOffset % 4)],
            seed.queryOffset, seed.subjectOffset);
	#endif

//	printf("Start=%d,%d End=%d,%d\n", beforeDpResults.best.queryOffset, beforeDpResults.best.subjectOffset,
//                                      afterDpResults.best.queryOffset, afterDpResults.best.subjectOffset);

    // Determine score by combining score from the two traces, and the match score at
    // the seed position
    return beforeDpResults.bestScore + afterDpResults.bestScore +
           choppedPSSMatrix.matrix[0][encoding_extractBase(choppedSubject[0],seed.subjectOffset % 4)];
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults nuGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                              struct coordinate seed, int4 dropoff)
{
    int2 **queryPosition, **bestQueryPosition;
    int2 **rowDropoff, **columnDropoff;
    unsigned char* subjectPosition, *bestSubjectPosition, subjectChar;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 queryDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char rightOfDropoff;
    struct dpResults dpResults;
    int4 bytePosition;

    if (seed.queryOffset == 0 || seed.subjectOffset == 0)
    {
        dpResults.best.queryOffset = 0;
        dpResults.best.subjectOffset = 0;
        dpResults.bestScore = bestScore;
        dpResults.traceback = NULL;
        return dpResults;
    }

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (seed.queryOffset >= nuGappedScoring_rowSizes)
    {
        // Free existing rows
        free(nuGappedScoring_bestRow);
        free(nuGappedScoring_insertQrow);
        // Set size to double current needed length
        nuGappedScoring_rowSizes = (seed.queryOffset) * 2;
        // Malloc new rows
        nuGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * nuGappedScoring_rowSizes);
        nuGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * nuGappedScoring_rowSizes);
    }

    // Convert subject offset to point to bytepacked subject
	bytePosition = (seed.subjectOffset - 1) % 4;
    bestSubjectPosition = subjectPosition = subject + ((seed.subjectOffset - 1) / 4);
    bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

    // Initialize row pointers
    rowOffset = (queryPosition - PSSMatrix.matrix);
    bestRow = nuGappedScoring_bestRow + rowOffset;
    insertQrow = nuGappedScoring_insertQrow + rowOffset;

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.matrix;
    columnDropoff = PSSMatrix.matrix + seed.queryOffset;

    // -----FIRST ROW-----
	subjectChar = encoding_extractBase(*subjectPosition, bytePosition);

    // -----FIRST CELL-----
    // Set M value for bottom-right cell
    match = (*queryPosition)[subjectChar];

    // M must be the best
    *bestRow = match;

    // Only gap opens possible
    *insertQrow = insertS = match - parameters_openGap;

    // If this is the best-yet scoring cell
    if (match > bestScore)
    {
        // Update best start cell data
        bestScore = match;
        bestQueryPosition = queryPosition;
        bestSubjectPosition = subjectPosition;
    }

    queryDistance = 0;
    queryPosition--; bestRow--; insertQrow--;

    // ----- REMAINING CELLS -----
    // For each remaining column in the bottom row, scanning from right-to-left
    while (queryPosition >= PSSMatrix.matrix)
    {
        // Set value for M
        match = (*queryPosition)[subjectChar]
              - parameters_openGap - queryDistance * parameters_extendGap;

        // Determine the best of M and Iy
        if (match > insertS)
        {
            *bestRow = match;

            // Calculate new Iy
            insertS = maximum(match - parameters_openGap,
                              insertS - parameters_extendGap);
        }
        else
        {
            *bestRow = insertS;

            // Since M <= Iy, new Iy must derive from Iy
            insertS -= parameters_extendGap;
        }

        // Set DUMMY Ix value, which should never be used
        *insertQrow = constants_gappedExtensionDummyValue;

        // If this is the best-yet scoring cell
        if (match > bestScore)
        {
            // Update best start cell data
            bestScore = match;
            bestQueryPosition = queryPosition;
            bestSubjectPosition = subjectPosition;
        }

        // If score at current cell is below dropoff
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = queryPosition;
            // And stop processing row
            break;
        }

        queryPosition--; bestRow--; insertQrow--;
        queryDistance++;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        nuGappedScoring_printBeforeRow(nuGappedScoring_bestRow, PSSMatrix.matrix, rowDropoff, columnDropoff);
    #endif

    // -----REMAINING ROWS-----
    while (rowDropoff < columnDropoff)
    {
    	// Move to next subject characters
        if (bytePosition)
        {
        	// Next char in current byte
        	bytePosition--;
        }
        else
        {
        	// Involves moving to next byte
        	bytePosition = 3;
            subjectPosition--;
            if (subjectPosition < subject)
            	break;
        }

		// Extract the subject characters
       	subjectChar = encoding_extractBase(*subjectPosition, bytePosition);

//    	printf("[%d/%d]", subjectPosition - subject, bytePosition);

        queryPosition = columnDropoff - 1;

        // Reset row pointers to start of rows
        rowOffset = (queryPosition - PSSMatrix.matrix);
        bestRow = nuGappedScoring_bestRow + rowOffset;
        insertQrow = nuGappedScoring_insertQrow + rowOffset;

        // -----FAR RIGHT CELL-----
        // Record some old values
        previousOldBest = *bestRow;

        // Ix is the best
        *bestRow = *insertQrow;

        // Calculate new Ix value
        *insertQrow -= parameters_extendGap;

        // Set DUMMY value for Iy, which should never be used
        insertS = constants_gappedExtensionDummyValue;

        // If score at current cell is below dropoff
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            columnDropoff = queryPosition;
            rightOfDropoff = 1;
        }
        else
        {
            // We are left of the column dropoff for this row
            rightOfDropoff = 0;
        }

        queryPosition--; bestRow--; insertQrow--;

        // -----CELLS RIGHT OF ROW DROPOFF-----
start1:
		// Loop 1 when insertS has no value
        while (queryPosition >= rowDropoff)
        {
            // Calculate new M value
            oldBest = *bestRow;
            match = (*queryPosition)[subjectChar] + previousOldBest;
            previousOldBest = oldBest;

            // Determine the best of M and Ix
            if (match > *insertQrow)
            {
                // Match is largest
                *bestRow = match;

                // Calculate new Ix
                *insertQrow = maximum(match - parameters_openGap,
                                      *insertQrow - parameters_extendGap);

                // Calculate new Iy
                insertS = maximum(match - parameters_openGap,
                                  insertS - parameters_extendGap);

                // If this is the best-yet scoring cell
                if (match > bestScore)
                {
                    // Update best start cell data
                    bestScore = match;
                    bestQueryPosition = queryPosition;
                    bestSubjectPosition = subjectPosition;
                }

                // If score at current cell (and cells to its right) are below dropoff
                if (rightOfDropoff)
                {
                    if (bestScore > *bestRow + dropoff)
                    {
                        // Record dropoff position
                        columnDropoff = queryPosition;
                    }
                    else
                    {
                        // We are left of the column dropoff for this row
                        rightOfDropoff = 0;
                    }
                }
                queryPosition--; bestRow--; insertQrow--;

                // InsertS now has a value
                break;
            }
            else
            {
                // insertQ is largest
                *bestRow = *insertQrow;

                // Calculate new Ix
                *insertQrow -= parameters_extendGap;
            }

            // If score at current cell (and cells to its right) are below dropoff
            if (rightOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = queryPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    rightOfDropoff = 0;
                }
            }

            queryPosition--; bestRow--; insertQrow--;
        }

        // Loop2 whilst insertS does have a value
        while (queryPosition >= rowDropoff)
        {
            // Calculate new M value
            oldBest = *bestRow;
            match = (*queryPosition)[subjectChar] + previousOldBest;
            previousOldBest = oldBest;

            // Determine the best of M, Ix and Iy
            if (match > insertS)
            {
            	if (match > *insertQrow)
                {
                	// Match is largest
                    *bestRow = match;

                    // Calculate new Ix
                    *insertQrow = maximum(match - parameters_openGap,
                                          *insertQrow - parameters_extendGap);

                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }
                }
                else
                {
                	// insertQ is largest
                    *bestRow = *insertQrow;

                    // Calculate new Ix
                    *insertQrow -= parameters_extendGap;

	                insertS = constants_gappedExtensionDummyValue;
                    break;
                }
            }
            else
            {
            	if (insertS > *insertQrow)
                {
                	// insertS is largest
                    *bestRow = insertS;

                    // Dummy Ix
                    *insertQrow = constants_gappedExtensionDummyValue;

                	// Calculate new Iy
	                insertS -= parameters_extendGap;

                }
                else
                {
                	// insertQ is largest
                    *bestRow = *insertQrow;

                    // Calculate new Ix
                    *insertQrow -= parameters_extendGap;

	                insertS = constants_gappedExtensionDummyValue;
                    break;
                }
            }

            // If score at current cell (and cells to its right) are below dropoff
            if (rightOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = queryPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    rightOfDropoff = 0;
                }
            }

            queryPosition--; bestRow--; insertQrow--;
        }

        if (queryPosition >= rowDropoff)
        {
            // If score at current cell (and cells to its right) are below dropoff
            if (rightOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = queryPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    rightOfDropoff = 0;
                }
            }

            queryPosition--; bestRow--; insertQrow--;
        	goto start1;
		}

        // -----CELLS LEFT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow + 1) + dropoff))
        {
            while (queryPosition >= PSSMatrix.matrix)
            {
                // Set value for Iy and best
                *bestRow = insertS;
                insertS = insertS - parameters_extendGap;

                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Stop processing row
                    queryPosition--;
                    break;
                }

                queryPosition--; bestRow--; insertQrow--;
            }
        }

        // Record dropoff position
        rowDropoff = queryPosition + 1;

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
            nuGappedScoring_printBeforeRow(nuGappedScoring_bestRow, PSSMatrix.matrix, rowDropoff, columnDropoff);
        #endif
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

// Perform dynamic programming to explore possible END points and alignments that start at
// the given seed and find the best score
struct dpResults nuGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                             int4 dropoff, int4 subjectLength)
{
    int2 **queryPosition, **bestQueryPosition, **queryEnd;
    int2 **rowDropoff, **columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd, subjectChar;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 queryDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char leftOfDropoff;
    int4 queryLength;
    struct dpResults dpResults;
    int4 bytePosition, subjectEndBytePosition;

    queryLength = PSSMatrix.length;
    queryEnd = PSSMatrix.matrix + queryLength;

    subjectEnd = subject + (subjectLength / 4);
    subjectEndBytePosition = subjectLength % 4;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (queryLength >= nuGappedScoring_rowSizes)
    {
        // Free existing rows
        free(nuGappedScoring_bestRow);
        free(nuGappedScoring_insertQrow);
        // Set size to double current needed length
        nuGappedScoring_rowSizes = queryLength * 2;
        // Malloc new rows
        nuGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * nuGappedScoring_rowSizes);
        nuGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * nuGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject;
	bytePosition = 1;
    bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

    // Initialize rows
    bestRow = nuGappedScoring_bestRow + 1;
    insertQrow = nuGappedScoring_insertQrow + 1;

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.matrix + queryLength - 1;
    columnDropoff = PSSMatrix.matrix;

    // -----FIRST ROW-----
	subjectChar = encoding_extractBase(*subjectPosition, bytePosition);

    // -----FIRST CELL-----
    // Set M value for top-left cell
    match = (*queryPosition)[subjectChar];

    // M must be the best
    *bestRow = match;

    // Only gap opens possible
    *insertQrow = insertS = match - parameters_openGap;

    // If this is the best-yet scoring cell
    if (match > bestScore)
    {
        // Update best start cell data
        bestScore = match;
        bestQueryPosition = queryPosition;
        bestSubjectPosition = subjectPosition;
    }

    queryDistance = 0;
    queryPosition++; bestRow++; insertQrow++;

    // ----- REMAINING CELLS -----
    // For each remaining columns in the top row, scanning from left-to-right
    while (queryPosition < queryEnd)
    {
        // Set value for M
        match = (*queryPosition)[subjectChar]
              - parameters_openGap - queryDistance * parameters_extendGap;

        // Determine the best of M and Iy
        if (match > insertS)
        {
            *bestRow = match;

            // Calculate new Iy
            insertS = maximum(match - parameters_openGap,
                              insertS - parameters_extendGap);
        }
        else
        {
            *bestRow = insertS;

            // Since M <= Iy, new Iy must derive from Iy
            insertS -= parameters_extendGap;
        }

        // Set DUMMY Ix value, which should never be used
        *insertQrow = constants_gappedExtensionDummyValue;

        // If this is the best-yet scoring cell
        if (match > bestScore)
        {
            // Update best start cell data
            bestScore = match;
            bestQueryPosition = queryPosition;
            bestSubjectPosition = subjectPosition;
        }

        // If score at current cell is below dropoff
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = queryPosition;
            // And stop processing row
            break;
        }

        queryPosition++; bestRow++; insertQrow++;
        queryDistance++;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
    	nuGappedScoring_printAfterRow(nuGappedScoring_bestRow + 1, PSSMatrix.matrix, rowDropoff, columnDropoff);
    #endif

    // Move to next subject character
    if (bytePosition == 3)
    {
    	bytePosition = 0;
	    subjectPosition++;
    }
    else
    {
    	bytePosition++;
    }

    // -----REMAINING ROWS-----
    while (rowDropoff > columnDropoff)
    {
    	// Stop at end of subject
    	if (subjectPosition == subjectEnd && bytePosition == subjectEndBytePosition)
        	break;

		subjectChar = encoding_extractBase(*subjectPosition, bytePosition);

        queryPosition = columnDropoff + 1;

        // Reset rows
        rowOffset = (queryPosition - PSSMatrix.matrix);
        bestRow = nuGappedScoring_bestRow + rowOffset;
        insertQrow = nuGappedScoring_insertQrow + rowOffset;

        // -----FAR LEFT CELL-----
        // Record some old values
        previousOldBest = *bestRow;

        // Ix is the best
        *bestRow = *insertQrow;

        // Calculate new Ix value
        *insertQrow -= parameters_extendGap;

        // Set DUMMY value for Iy, which should never be used
        insertS = constants_gappedExtensionDummyValue;

        // If score at current cell is below dropoff
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            columnDropoff = queryPosition;
            leftOfDropoff = 1;
        }
        else
        {
            // We are left of the column dropoff for this row
            leftOfDropoff = 0;
        }

        queryPosition++; bestRow++; insertQrow++;

        // -----CELLS LEFT OF ROW DROPOFF-----
start2:
		// Loop 1 when insertS has no value
        while (queryPosition <= rowDropoff)
        {
            // Calculate new M value
            oldBest = *bestRow;
            match = (*queryPosition)[subjectChar] + previousOldBest;
            previousOldBest = oldBest;

            // Determine the best of M and Ix
            if (match > *insertQrow)
            {
                // Match is largest
                *bestRow = match;

                // Calculate new Ix
                *insertQrow = maximum(match - parameters_openGap,
                                      *insertQrow - parameters_extendGap);

                // Calculate new Iy
                insertS = maximum(match - parameters_openGap,
                                  insertS - parameters_extendGap);

                // If this is the best-yet scoring cell
                if (match > bestScore)
                {
                    // Update best start cell data
                    bestScore = match;
                    bestQueryPosition = queryPosition;
                    bestSubjectPosition = subjectPosition;
                }

                // If score at current cell (and cells to its right) are below dropoff
                if (leftOfDropoff)
                {
                    if (bestScore > *bestRow + dropoff)
                    {
                        // Record dropoff position
                        columnDropoff = queryPosition;
                    }
                    else
                    {
                        // We are left of the column dropoff for this row
                        leftOfDropoff = 0;
                    }
                }
                queryPosition++; bestRow++; insertQrow++;

                // InsertS now has a value
                break;
            }
            else
            {
                // insertQ is largest
                *bestRow = *insertQrow;

                // Calculate new Ix
                *insertQrow -= parameters_extendGap;
            }

            // If score at current cell (and cells to its right) are below dropoff
            if (leftOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = queryPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    leftOfDropoff = 0;
                }
            }

            queryPosition++; bestRow++; insertQrow++;
        }

        // Loop2 whilst insertS does have a value
        while (queryPosition <= rowDropoff)
        {
            // Calculate new M value
            oldBest = *bestRow;
            match = (*queryPosition)[subjectChar] + previousOldBest;
            previousOldBest = oldBest;

            // Determine the best of M, Ix and Iy
            if (match > insertS)
            {
            	if (match > *insertQrow)
                {
                	// Match is largest
                    *bestRow = match;

                    // Calculate new Ix
                    *insertQrow = maximum(match - parameters_openGap,
                                          *insertQrow - parameters_extendGap);

                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }
                }
                else
                {
                	// insertQ is largest
                    *bestRow = *insertQrow;

                    // Calculate new Ix
                    *insertQrow -= parameters_extendGap;

	                insertS = constants_gappedExtensionDummyValue;
                    break;
                }
            }
            else
            {
            	if (insertS > *insertQrow)
                {
                	// insertS is largest
                    *bestRow = insertS;

                    // Dummy Ix
                    *insertQrow = constants_gappedExtensionDummyValue;

                	// Calculate new Iy
	                insertS -= parameters_extendGap;

                }
                else
                {
                	// insertQ is largest
                    *bestRow = *insertQrow;

                    // Calculate new Ix
                    *insertQrow -= parameters_extendGap;

	                insertS = constants_gappedExtensionDummyValue;
                    break;
                }
            }

            // If score at current cell (and cells to its right) are below dropoff
            if (leftOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = queryPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    leftOfDropoff = 0;
                }
            }

            queryPosition++; bestRow++; insertQrow++;
        }

        if (queryPosition <= rowDropoff)
        {
            // If score at current cell (and cells to its right) are below dropoff
            if (leftOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = queryPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    leftOfDropoff = 0;
                }
            }

            queryPosition++; bestRow++; insertQrow++;
        	goto start2;
		}

        // -----CELLS RIGHT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow - 1) + dropoff))
        {
            while (queryPosition < queryEnd)
            {
                // Set value for Iy and best
                *bestRow = insertS;
                insertS = insertS - parameters_extendGap;

                // Set DUMMY value for Ix, which should never be used
                *insertQrow = constants_gappedExtensionDummyValue;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // And stop processing row
                    queryPosition++;
                    break;
                }

                queryPosition++; bestRow++; insertQrow++;
            }
        }

        // Record dropoff position
        rowDropoff = queryPosition - 1;

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

        // Move to next subject character
        if (bytePosition == 3)
        {
            bytePosition = 0;
            subjectPosition++;
        }
        else
        {
            bytePosition++;
        }

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
            nuGappedScoring_printAfterRow(nuGappedScoring_bestRow + 1, PSSMatrix.matrix, rowDropoff, columnDropoff);
        #endif
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

void nuGappedScoring_free()
{
    free(nuGappedScoring_bestRow);
    free(nuGappedScoring_insertQrow);
}

// Debugging routine
void nuGappedScoring_printBeforeRow(int4* row, int2** query, int2** rowDropoff, int2** columnDropoff)
{
	int2** queryPosition = query;
    int4 *oldRow, *best;

	while (queryPosition < rowDropoff)
	{
		printf("     ");
		queryPosition++;
		row++;
	}

    best = oldRow = row;

	while (queryPosition < columnDropoff)
	{
    	if (*row > *best)
        	best = row;
		row++;
		queryPosition++;
	}

    row = oldRow; queryPosition = rowDropoff;

    while (queryPosition < columnDropoff)
	{
    	if (row == best)
			printf("%4d*", *row);
        else
			printf("%4d ", *row);
		row++;
		queryPosition++;
	}

	printf("\n");
}

// Debugging routine
void nuGappedScoring_printAfterRow(int4* row, int2** query, int2** rowDropoff, int2** columnDropoff)
{
	int2** queryPosition = query;
    int4 *oldRow, *best;

	while (queryPosition < columnDropoff)
	{
		printf("     ");
		queryPosition++;
		row++;
	}

    best = oldRow = row;

	while (queryPosition < rowDropoff)
	{
    	if (*row > *best)
        	best = row;
		row++;
		queryPosition++;
	}

    row = oldRow; queryPosition = columnDropoff;

    while (queryPosition < rowDropoff)
	{
    	if (row == best)
			printf("%4d*", *row);
        else
			printf("%4d ", *row);
		row++;
		queryPosition++;
	}

	printf("\n");
}

