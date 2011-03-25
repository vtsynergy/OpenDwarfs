// semiGappedScoring.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform semi-gapped alignment WITHOUT restricted insertion
// recording the score and alignment end-points only

#include "blast.h"

int4* oldSemiGappedScoring_bestRow = NULL;
int4* oldSemiGappedScoring_insertQrow = NULL;
int4 oldSemiGappedScoring_rowSizes = 0;

#define maximum(a,b) ((a > b) ? a : b)

// Prototypes
struct dpResults oldSemiGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                               struct coordinate seed, int4 dropoff);
struct dpResults oldSemiGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                              int4 dropoff, int4 subjectLength);

// Perform semi-gapped alignment without restricted insertion
int4 oldSemiGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
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

    beforeDpResults = oldSemiGappedScoring_dpBeforeSeed(subject, PSSMatrix,
                    seed, dropoff + parameters_semiGappedDropoffIncrease);

    // Chop the start off the query and subject so they begin at the seed
    choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
    choppedSubject = subject + seed.subjectOffset;
    choppedSubjectSize = subjectSize - seed.subjectOffset;

    // Perform dynamic programming for points after the seed
    afterDpResults = oldSemiGappedScoring_dpAfterSeed(choppedSubject, choppedPSSMatrix,
                    dropoff + parameters_semiGappedDropoffIncrease, choppedSubjectSize);

    // Re-adjust result change due to chopping subject/query and strand adjustment
    afterDpResults.best.queryOffset += seed.queryOffset + strandOffset;
    afterDpResults.best.subjectOffset += seed.subjectOffset;
    beforeDpResults.best.queryOffset += strandOffset;

    // Associate best scoring start and end points with the ungapped extension
    ungappedExtension->start = beforeDpResults.best;
    ungappedExtension->end = afterDpResults.best;

//    if (dloc == 88197331)
//        printf("semiGapped[%d,%d,%d]\n", beforeDpResults.bestScore, afterDpResults.bestScore,
//        choppedPSSMatrix.matrix[0][choppedSubject[0]]);

    // Determine score by combining score from the two traces, and the match score at
    // the seed position
    return beforeDpResults.bestScore + afterDpResults.bestScore +
           choppedPSSMatrix.matrix[0][choppedSubject[0]];
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults oldSemiGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                 struct coordinate seed, int4 dropoff)
{
    int2 **queryPosition, **bestQueryPosition;
    int2* matrixColumn;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *startSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 subjectDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char rightOfDropoff;
    int4 queryCount, subjectCount;
    struct dpResults dpResults;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (seed.subjectOffset >= oldSemiGappedScoring_rowSizes)
    {
        // Free existing rows
        free(oldSemiGappedScoring_bestRow);
        free(oldSemiGappedScoring_insertQrow);
        // Set size to double current needed length
        oldSemiGappedScoring_rowSizes = (seed.subjectOffset) * 2;
        // Malloc new rows
        oldSemiGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * oldSemiGappedScoring_rowSizes);
        oldSemiGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * oldSemiGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = startSubjectPosition = subject + seed.subjectOffset - 1;
    bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

    // Initialize row pointers
    rowOffset = (subjectPosition - subject);
    bestRow = oldSemiGappedScoring_bestRow + rowOffset;
    insertQrow = oldSemiGappedScoring_insertQrow + rowOffset;

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
    *insertQrow = insertS = match - parameters_semiGappedOpenGap;

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
              - parameters_semiGappedOpenGap - subjectDistance * parameters_semiGappedExtendGap;

        // Determine the best of M and Iy
        if (match > insertS)
        {
            *bestRow = match;

            // Calculate new Iy
            insertS = maximum(match - parameters_semiGappedOpenGap,
                              insertS - parameters_semiGappedExtendGap);
        }
        else
        {
            *bestRow = insertS;

            // Since M <= Iy, new Iy must derive from Iy
            insertS -= parameters_semiGappedExtendGap;
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

//    if (dloc == 746829265)
//    print(oldSemiGappedScoring_bestRow, subject, rowDropoff, columnDropoff);

    // Start queryCount at N. Only allow insertS for every Nth row when queryCount
    // reaches 0
    queryCount = parameters_semiGappedExtensionN;

    // -----REMAINING ROWS-----
    while (queryPosition > PSSMatrix.matrix && rowDropoff < columnDropoff)
    {
        queryPosition--;
        queryCount--;
        subjectPosition = columnDropoff - 1;

        // Determine subjectCount for initial subjectPosition. Is used to only allow
        // insertQ when subjectOffset % parameters_semiGappedExtensionN == 0
        subjectCount = (int4)(startSubjectPosition - subjectPosition) % parameters_semiGappedExtensionN;
        if (subjectCount)
            subjectCount = parameters_semiGappedExtensionN - subjectCount;

        // Reset row pointers to start of rows
        rowOffset = (subjectPosition - subject);
        bestRow = oldSemiGappedScoring_bestRow + rowOffset;
        insertQrow = oldSemiGappedScoring_insertQrow + rowOffset;

        // Using next column of query matrix
        matrixColumn = *queryPosition;

        // ************ All rows we are not allowing insertS
        if (queryCount)
        {
            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                rightOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // -----FAR RIGHT CELL-----
                // Record some old values
                previousOldBest = *bestRow;

                // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

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

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition--; bestRow--; insertQrow--; subjectCount--;

            // -----CELLS RIGHT OF ROW DROPOFF-----
            while (subjectPosition >= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Calculate new M value, which is also the best
                    oldBest = *bestRow;
                    match = *bestRow = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }
                }
                // We are allowing insertQ this column
                else
                {
                    // Calculate new M value
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Ix
					if (match > *insertQrow)
                    {
						*bestRow = match;

                        // Calculate new Ix
                        *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                              *insertQrow - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = *insertQrow;

                        // Since M <= Ix, new Ix must derive from Ix
						*insertQrow -= parameters_semiGappedExtendGap;
                    }

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition--; bestRow--; insertQrow--; subjectCount--;
            }

            // -----SINGLE CELL LEFT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition >= subject))
            {
                // Set value for best
                *bestRow = match = previousOldBest + matrixColumn[*subjectPosition];

                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                if (match + dropoff >= bestScore)
                {
                    // Record dropoff position
                    rowDropoff = subjectPosition;
                }
            }
        }

        // ************ Every Nth row we allow insertS
        else
        {
            // -----FAR RIGHT CELL-----

            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                rightOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // Record some old values
                previousOldBest = *bestRow;

                // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

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

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition--; bestRow--; insertQrow--; subjectCount--;

            // -----CELLS RIGHT OF ROW DROPOFF-----
            while (subjectPosition >= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Iy
					if (match > insertS)
                    {
						*bestRow = match;

                        // Calculate new Iy
                        insertS = maximum(match - parameters_semiGappedOpenGap,
                                          insertS - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = insertS;

                        // Since M <= Iy, new Iy must derive from Iy
						insertS -= parameters_semiGappedExtendGap;
                    }

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
                            columnDropoff = subjectPosition;
                        }
                        else
                        {
                            // We are left of the column dropoff for this row
                            rightOfDropoff = 0;
                        }
                    }
                }
                // ** We are allowing insertQ this column
                else
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
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
                            *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                                  *insertQrow - parameters_semiGappedExtendGap);

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);

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
                            *insertQrow -= parameters_semiGappedExtendGap;

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);
                        }
                    }
                    else
                    {
                        if (insertS > *insertQrow)
                        {
                            // insertS is largest
                            *bestRow = insertS;

                            // Calculate new Ix
                            *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                                  *insertQrow - parameters_semiGappedExtendGap);

                            // Calculate new Iy
                            insertS -= parameters_semiGappedExtendGap;
                        }
                        else
                        {
                            // insertQ is largest
                            *bestRow = *insertQrow;

                            // Calculate new Ix
                            *insertQrow -= parameters_semiGappedExtendGap;

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);
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

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition--; bestRow--; insertQrow--; subjectCount--;
            }

            // -----SINGLE CELL LEFT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition >= subject))
            {
                // Calculate match value
                match = previousOldBest + matrixColumn[*subjectPosition];

                // Set value for best
                *bestRow = maximum(match, insertS);

                // Calculate new Iy
                insertS = maximum(match - parameters_semiGappedOpenGap,
                                  insertS - parameters_semiGappedExtendGap);
 
                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                subjectPosition--; bestRow--; insertQrow--;
            }

            // -----CELLS LEFT OF ROW DROPOFF -----
            if (!(bestScore > *(bestRow + 1) + dropoff))
            {
                while (subjectPosition >= subject)
                {
                    // Set value for Iy and best
                    *bestRow = insertS;
                    insertS = insertS - parameters_semiGappedExtendGap;

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
                }
            }

            // Record dropoff position
            rowDropoff = subjectPosition + 1;

            // Clear insertS for next row
            insertS = constants_gappedExtensionDummyValue;

            // Reset queryCount
            queryCount = parameters_semiGappedExtensionN;
        }
//        if (dloc == 746829265)
//                print(oldSemiGappedScoring_bestRow, subject, rowDropoff, columnDropoff);
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

// Perform dynamic programming to explore possible END points and alignments that start at
// the given seed and find the best score
struct dpResults oldSemiGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                              int4 dropoff, int4 subjectLength)
{
    int2 **queryPosition, **bestQueryPosition, **queryEnd;
    int2* matrixColumn;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd, *startSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 subjectDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char leftOfDropoff;
    int4 queryLength;
    int4 queryCount, subjectCount;
    struct dpResults dpResults;

    queryLength = PSSMatrix.length;
    subjectEnd = subject + subjectLength;
    queryEnd = PSSMatrix.matrix + queryLength;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (subjectLength >= oldSemiGappedScoring_rowSizes)
    {
        // Free existing rows
        free(oldSemiGappedScoring_bestRow);
        free(oldSemiGappedScoring_insertQrow);
        // Set size to double current needed length
        oldSemiGappedScoring_rowSizes = subjectLength * 2;
        // Malloc new rows
        oldSemiGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * oldSemiGappedScoring_rowSizes);
        oldSemiGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * oldSemiGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = startSubjectPosition = subject + 1;
    bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

    // Initialize rows
    bestRow = oldSemiGappedScoring_bestRow + 1;
    insertQrow = oldSemiGappedScoring_insertQrow + 1;

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
    *insertQrow = insertS = match - parameters_semiGappedOpenGap;

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
              - parameters_semiGappedOpenGap - subjectDistance * parameters_semiGappedExtendGap;

        // Determine the best of M and Iy
        if (match > insertS)
        {
            *bestRow = match;

            // Calculate new Iy
            insertS = maximum(match - parameters_semiGappedOpenGap,
                              insertS - parameters_semiGappedExtendGap);
        }
        else
        {
            *bestRow = insertS;

            // Since M <= Iy, new Iy must derive from Iy
            insertS -= parameters_semiGappedExtendGap;
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

//    if (dloc == 88197331)
//        print2(oldSemiGappedScoring_bestRow + 1, subject, rowDropoff, columnDropoff);

    // Start queryCount at N. Only allow insertS for every Nth row when queryCount
    // reaches 0
    queryCount = parameters_semiGappedExtensionN;

    queryPosition++; queryCount--;

    // -----REMAINING ROWS-----
    while (queryPosition < queryEnd && rowDropoff > columnDropoff)
    {
        subjectPosition = columnDropoff + 1;

        // Determine subjectCount for initial subjectPosition. Is used to only allow
        // insertQ when subjectOffset % parameters_semiGappedExtensionN == 0
        subjectCount = ((int4)(subjectPosition - startSubjectPosition) % parameters_semiGappedExtensionN);
        if (subjectCount)
            subjectCount = parameters_semiGappedExtensionN - subjectCount;

        // Reset rows
        rowOffset = (subjectPosition - subject);
        bestRow = oldSemiGappedScoring_bestRow + rowOffset;
        insertQrow = oldSemiGappedScoring_insertQrow + rowOffset;

        // Using next column of query matrix
        matrixColumn = *queryPosition;


        // ************ All rows we are not allowing insertS
        if (queryCount)
        {
            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                leftOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // -----FAR LEFT CELL-----
                // Record some old values
                previousOldBest = *bestRow;

                // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = subjectPosition;
                    leftOfDropoff = 1;
                }
                else
                {
                    // We are right of the column dropoff for this row
                    leftOfDropoff = 0;
                }

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition++; bestRow++; insertQrow++; subjectCount--;

            // -----CELLS LEFT OF ROW DROPOFF-----
            while (subjectPosition <= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Calculate new M value, which is also the best
                    oldBest = *bestRow;
                    match = *bestRow = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }
                }
                // We are allowing insertQ this column
                else
                {
                    // Calculate new M value
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Ix
					if (match > *insertQrow)
                    {
						*bestRow = match;

                        // Calculate new Ix
                        *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                              *insertQrow - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = *insertQrow;

                        // Since M <= Ix, new Ix must derive from Ix
						*insertQrow -= parameters_semiGappedExtendGap;
                    }

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition++; bestRow++; insertQrow++; subjectCount--;
            }

            // -----SINGLE CELL RIGHT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition < subjectEnd))
            {
                // Set value for best
                *bestRow = match = previousOldBest + matrixColumn[*subjectPosition];

                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                if (match + dropoff >= bestScore)
                {
                    // Record dropoff position
                    rowDropoff = subjectPosition;
                }
            }
        }

        // ************ Every Nth row we allow insertS
        else
        {
            // -----FAR LEFT CELL-----

            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                leftOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // Record some old values
                previousOldBest = *bestRow;

                 // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

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
                    // We are right of the column dropoff for this row
                    leftOfDropoff = 0;
                }

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition++; bestRow++; insertQrow++; subjectCount--;

            // -----CELLS LEFT OF ROW DROPOFF-----
            while (subjectPosition <= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Iy
					if (match > insertS)
                    {
						*bestRow = match;

                        // Calculate new Iy
                        insertS = maximum(match - parameters_semiGappedOpenGap,
                                          insertS - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = insertS;

                        // Since M <= Iy, new Iy must derive from Iy
						insertS -= parameters_semiGappedExtendGap;
                    }

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
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
                            // We are right of the column dropoff for this row
                            leftOfDropoff = 0;
                        }
                    }
                }
                // ** We are allowing insertQ this column
                else
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
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
                            *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                                  *insertQrow - parameters_semiGappedExtendGap);

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);

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
                            *insertQrow -= parameters_semiGappedExtendGap;

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);
                        }
                    }
                    else
                    {
                        if (insertS > *insertQrow)
                        {
                            // insertS is largest
                            *bestRow = insertS;

                            // Calculate new Ix
                            *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                                  *insertQrow - parameters_semiGappedExtendGap);

                            // Calculate new Iy
                            insertS -= parameters_semiGappedExtendGap;
                        }
                        else
                        {
                            // insertQ is largest
                            *bestRow = *insertQrow;

                            // Calculate new Ix
                            *insertQrow -= parameters_semiGappedExtendGap;

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);
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
                            // We are right of the column dropoff for this row
                            leftOfDropoff = 0;
                        }
                    }

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition++; bestRow++; insertQrow++; subjectCount--;
            }

            // -----SINGLE CELL RIGHT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition < subjectEnd))
            {
                // Calculate match value
                match = previousOldBest + matrixColumn[*subjectPosition];

                // Set value for best
                *bestRow = maximum(match, insertS);

                // Calculate new Iy
                insertS = maximum(match - parameters_semiGappedOpenGap,
                                  insertS - parameters_semiGappedExtendGap);

                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                subjectPosition++; bestRow++; insertQrow++;
            }

            // -----CELLS RIGHT OF ROW DROPOFF -----
            if (!(bestScore > *(bestRow - 1) + dropoff))
            {
                while (subjectPosition < subjectEnd)
                {
                    // Set value for Iy and best
                    *bestRow = insertS;
                    insertS = insertS - parameters_semiGappedExtendGap;

                    // Set DUMMY values for Ix
                    *insertQrow = constants_gappedExtensionDummyValue;

                    // If score at current cell is below dropoff
                    if (bestScore > *bestRow + dropoff)
                    {
                        // Stop processing row
                        subjectPosition++;
                        break;
                    }

                    subjectPosition++; bestRow++; insertQrow++;
                }
            }

            // Record dropoff position
            rowDropoff = subjectPosition - 1;

            // Clear insertS for next row
            insertS = constants_gappedExtensionDummyValue;

            // Reset queryCount
            queryCount = parameters_semiGappedExtensionN;
        }
//    if (dloc == 88197331)
//        print2(oldSemiGappedScoring_bestRow + 1, subject, rowDropoff, columnDropoff);

        queryPosition++; queryCount--;
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

void oldSemiGappedScoring_free()
{
    free(oldSemiGappedScoring_bestRow);
    free(oldSemiGappedScoring_insertQrow);
}
