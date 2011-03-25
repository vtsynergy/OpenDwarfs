// oldGappedScoring.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform traditional gapped alignment WITHOUT restricted insertion
// recording the score and alignment end-points only

#include "blast.h"

int4* oldGappedScoring_bestRow = NULL;
int4* oldGappedScoring_insertQrow = NULL;
int4 oldGappedScoring_rowSizes = 0;

#define maximum(a,b) ((a > b) ? a : b)

// Prototypes
struct dpResults oldGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                               struct coordinate seed, int4 dropoff);
struct dpResults oldGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                              int4 dropoff, int4 subjectLength);

int4 oldGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff)
{
    struct coordinate seed;
    unsigned char *choppedSubject;
    struct PSSMatrix choppedPSSMatrix;
    int4 choppedSubjectSize;
    struct dpResults beforeDpResults, afterDpResults;
    int4 strandOffset = 0;

    // Perform dynamic programming for points before the seed
    seed = ungappedExtension->seed;
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

//    ungappedExtension_print(ungappedExtension);
//    printf("dloc=%d\n", blast_dloc);
//    printf("seed=%d,%d\n", seed.queryOffset, seed.subjectOffset); fflush(stdout);

    beforeDpResults = oldGappedScoring_dpBeforeSeed(subject, PSSMatrix, seed, dropoff);

    // Chop the start off the query and subject so they begin at the seed
    choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
    choppedSubject = subject + seed.subjectOffset;
    choppedSubjectSize = subjectSize - seed.subjectOffset;

    // Perform dynamic programming for points after the seed
    afterDpResults = oldGappedScoring_dpAfterSeed(choppedSubject, choppedPSSMatrix,
                                               dropoff, choppedSubjectSize);

    // Re-adjust result change due to chopping subject/query and strand adjustment
    afterDpResults.best.queryOffset += seed.queryOffset + strandOffset;
    afterDpResults.best.subjectOffset += seed.subjectOffset;
    beforeDpResults.best.queryOffset += strandOffset;

    // Associate best scoring start and end points with the ungapped extension
    ungappedExtension->start = beforeDpResults.best;
    ungappedExtension->end = afterDpResults.best;

//    if (dloc == 1338183)
//        printf("norm[%d,%d,%d](%d)\n", beforeDpResults.bestScore, afterDpResults.bestScore,
//        choppedPSSMatrix.matrix[0][choppedSubject[0]], seed.queryOffset);

    // Determine score by combining score from the two traces, and the match score at
    // the seed position
    return beforeDpResults.bestScore + afterDpResults.bestScore +
           choppedPSSMatrix.matrix[0][choppedSubject[0]];
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults oldGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                            struct coordinate seed, int4 dropoff)
{
    int2 **queryPosition, **bestQueryPosition;
    int2* matrixColumn;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char* subjectPosition, *bestSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 subjectDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char rightOfDropoff;
    struct dpResults dpResults;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (seed.subjectOffset >= oldGappedScoring_rowSizes)
    {
        // Free existing rows
        free(oldGappedScoring_bestRow);
        free(oldGappedScoring_insertQrow);
        // Set size to double current needed length
        oldGappedScoring_rowSizes = (seed.subjectOffset) * 2;
        // Malloc new rows
        oldGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * oldGappedScoring_rowSizes);
        oldGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * oldGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject + seed.subjectOffset - 1;
    bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

    // Initialize row pointers
    rowOffset = (subjectPosition - subject);
    bestRow = oldGappedScoring_bestRow + rowOffset;
    insertQrow = oldGappedScoring_insertQrow + rowOffset;

    // Set initial row dropoff and column dropoff
    rowDropoff = subject;
    columnDropoff = subject + seed.subjectOffset;

    // Using first column of query matrix
    matrixColumn = *queryPosition;

    // -----FIRST ROW-----

    // -----FIRST CELL-----
    // Set M value for bottom-right cell
    match = matrixColumn[*subjectPosition];

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

    subjectDistance = 0;
    subjectPosition--; bestRow--; insertQrow--;

    // ----- REMAINING CELLS -----
    // For each remaining column in the bottom row, scanning from right-to-left
    while (subjectPosition >= subject)
    {
        // Set value for M
        match = matrixColumn[*subjectPosition]
              - parameters_openGap - subjectDistance * parameters_extendGap;

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
            rowDropoff = subjectPosition;
            // And stop processing row
            break;
        }

        subjectPosition--; bestRow--; insertQrow--;
        subjectDistance++;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

//    if (dloc == 19063576)
//    print(oldGappedScoring_bestRow, subject, rowDropoff, columnDropoff);

    // -----REMAINING ROWS-----
    while (queryPosition > PSSMatrix.matrix && rowDropoff < columnDropoff)
    {
        queryPosition--;
        subjectPosition = columnDropoff - 1;

        // Reset row pointers to start of rows
        rowOffset = (subjectPosition - subject);
        bestRow = oldGappedScoring_bestRow + rowOffset;
        insertQrow = oldGappedScoring_insertQrow + rowOffset;

        // Using next column of query matrix
        matrixColumn = *queryPosition;

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
            columnDropoff = subjectPosition;
            rightOfDropoff = 1;
        }
        else
        {
            // We are left of the column dropoff for this row
            rightOfDropoff = 0;
        }

        subjectPosition--; bestRow--; insertQrow--;

        // -----CELLS RIGHT OF ROW DROPOFF-----
        while (subjectPosition >= rowDropoff)
        {
            // Remember old M value (for cell below this one)
            oldBest = *bestRow;

            // Calculate new M value
            match = matrixColumn[*subjectPosition] + previousOldBest;
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

                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);
                }
            }
            else
            {
            	if (insertS > *insertQrow)
                {
                	// insertS is largest
                    *bestRow = insertS;

                    // Calculate new Ix
                    *insertQrow = maximum(match - parameters_openGap,
                                          *insertQrow - parameters_extendGap);

                	// Calculate new Iy
	                insertS -= parameters_extendGap;

                }
                else
                {
                	// insertQ is largest
                    *bestRow = *insertQrow;

                    // Calculate new Ix
                    *insertQrow -= parameters_extendGap;

                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);
                }
            }

            // If score at current cell (and cells to its right) are below dropoff
            if (rightOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = subjectPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    rightOfDropoff = 0;
                }
            }

            subjectPosition--; bestRow--; insertQrow--;
        }

        // -----CELLS LEFT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow + 1) + dropoff))
        {
            while (subjectPosition >= subject)
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
                    subjectPosition--;
                    break;
                }

                subjectPosition--; bestRow--; insertQrow--;
                subjectDistance++;
            }
        }

        // Record dropoff position
        rowDropoff = subjectPosition + 1;

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

//    if (dloc == 20877970)
//        print(oldGappedScoring_bestRow, subject, rowDropoff, columnDropoff);
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

// Perform dynamic programming to explore possible END points and alignments that start at
// the given seed and find the best score
struct dpResults oldGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                           int4 dropoff, int4 subjectLength)
{
    int2 **queryPosition, **bestQueryPosition, **queryEnd;
    int2* matrixColumn;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 subjectDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char leftOfDropoff;
    int4 queryLength;
    struct dpResults dpResults;

    queryLength = PSSMatrix.length;
    subjectEnd = subject + subjectLength;
    queryEnd = PSSMatrix.matrix + queryLength;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (subjectLength >= oldGappedScoring_rowSizes)
    {
        // Free existing rows
        free(oldGappedScoring_bestRow);
        free(oldGappedScoring_insertQrow);
        // Set size to double current needed length
        oldGappedScoring_rowSizes = subjectLength * 2;
        // Malloc new rows
        oldGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * oldGappedScoring_rowSizes);
        oldGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * oldGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject + 1;
    bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

    // Initialize rows
    bestRow = oldGappedScoring_bestRow + 1;
    insertQrow = oldGappedScoring_insertQrow + 1;

    // Set initial row dropoff and column dropoff
    rowDropoff = subject + subjectLength - 1;
    columnDropoff = subject;

    // -----FIRST ROW-----

    // Using first column of query matrix
    matrixColumn = *queryPosition;

    // -----FIRST CELL-----
    // Set M value for top-left cell
    match = matrixColumn[*subjectPosition];

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

    subjectDistance = 0;
    subjectPosition++; bestRow++; insertQrow++;

    // ----- REMAINING CELLS -----
    // For each remaining columns in the top row, scanning from left-to-right
    while (subjectPosition < subjectEnd)
    {
        // Set value for M
        match = matrixColumn[*subjectPosition]
              - parameters_openGap - subjectDistance * parameters_extendGap;

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
            rowDropoff = subjectPosition;
            // And stop processing row
            break;
        }

        subjectPosition++; bestRow++; insertQrow++;
        subjectDistance++;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

//    if (dloc == 94501334)
//    print2(oldGappedScoring_bestRow + 1, subject, rowDropoff, columnDropoff);

    queryPosition++;

    // -----REMAINING ROWS-----
    while (queryPosition < queryEnd && rowDropoff > columnDropoff)
    {
        subjectPosition = columnDropoff + 1;

        // Using next column of query matrix
        matrixColumn = *queryPosition;

        // Reset rows
        rowOffset = (subjectPosition - subject);
        bestRow = oldGappedScoring_bestRow + rowOffset;
        insertQrow = oldGappedScoring_insertQrow + rowOffset;

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
            columnDropoff = subjectPosition;
            leftOfDropoff = 1;
        }
        else
        {
            // We are left of the column dropoff for this row
            leftOfDropoff = 0;
        }

        subjectPosition++; bestRow++; insertQrow++;

        // -----CELLS LEFT OF ROW DROPOFF-----
        while (subjectPosition <= rowDropoff)
        {
            // Remember old Best value (for cell below this one)
            oldBest = *bestRow;

            // Calculate new M value
            match = matrixColumn[*subjectPosition] + previousOldBest;

            // Record some old values
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

                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);
                }
            }
            else
            {
            	if (insertS > *insertQrow)
                {
                	// insertS is largest
                    *bestRow = insertS;

                    // Calculate new Ix
                    *insertQrow = maximum(match - parameters_openGap,
                                          *insertQrow - parameters_extendGap);

                	// Calculate new Iy
	                insertS -= parameters_extendGap;

                }
                else
                {
                	// insertQ is largest
                    *bestRow = *insertQrow;

                    // Calculate new Ix
                    *insertQrow -= parameters_extendGap;

                    // Calculate new Iy
                    insertS = maximum(match - parameters_openGap,
                                      insertS - parameters_extendGap);
                }
            }

            // If score at current cell (and cells to its left) are below dropoff
            if (leftOfDropoff)
            {
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = subjectPosition;
                }
                else
                {
                    // We are left of the column dropoff for this row
                    leftOfDropoff = 0;
                }
            }

            subjectPosition++; bestRow++; insertQrow++;
        }

        // -----CELLS RIGHT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow - 1) + dropoff))
        {
            while (subjectPosition < subjectEnd)
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
                    subjectPosition++;
                    break;
                }

                subjectPosition++; bestRow++; insertQrow++;
                subjectDistance++;
            }
        }

        // Record dropoff position
        rowDropoff = subjectPosition - 1;

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

        queryPosition++;

//    if (dloc == 19063576)
//        print2(oldGappedScoring_bestRow + 1, subject, rowDropoff, columnDropoff);
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

void oldGappedScoring_free()
{
    free(oldGappedScoring_bestRow);
    free(oldGappedScoring_insertQrow);
}
