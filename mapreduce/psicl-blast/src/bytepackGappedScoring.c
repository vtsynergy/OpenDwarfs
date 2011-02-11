// bytepackGappedScoring.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform gapped alignment using byte-packed alignment technique

#include "blast.h"
#include "fasterBytepackGappedScoring.c"

int4* bytepackGappedScoring_bestRows[2] = {NULL,NULL};
int4* bytepackGappedScoring_insertQrow = NULL;
int4 bytepackGappedScoring_rowSizes = 0;

#define maximum(a,b) ((a > b) ? a : b)

// Prototypes
struct dpResults bytepackGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                    struct coordinate seed, int4 dropoff);
struct dpResults bytepackGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                   int4 dropoff, int4 subjectLength);

int4 bytepackGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                                int4 subjectSize, unsigned char* packedSubject, int4 dropoff)
{
    struct coordinate seed;
    unsigned char *choppedSubject;
    struct PSSMatrix choppedPSSMatrix;
    int4 choppedSubjectSize, score, packedSubjectSize;
    struct dpResults beforeDpResults, afterDpResults;
    int4 strandOffset = 0;

    // Increment counter number of HSP's semi-gapped
	blast_numSemiGapped++;

    dropoff -= parameters_bytepackDropoffDecrease;

    // Calculate packed subject size
    packedSubjectSize = (subjectSize + 3) / 4;

    seed = ungappedExtension->seed;
    if (seed.queryOffset > PSSMatrix.strandLength)
    {
	    // If query position is in the second strand, remove first strand from PSSM
        strandOffset = PSSMatrix.strandLength;
		seed.queryOffset -= PSSMatrix.strandLength;
		PSSMatrix = PSSMatrix_chop(PSSMatrix, PSSMatrix.strandLength);

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
			printf("Second stand qLength=%d seed=%d,%d\n", PSSMatrix.length,
                   seed.queryOffset, seed.subjectOffset);
        #endif
    }
    else
    {
    	// Otherwise remove second strand
    	PSSMatrix.length = PSSMatrix.strandLength;
    }

    // Adjust to point to start of byte
    seed.queryOffset -= 2;
    seed.subjectOffset -= 2;

//    printf("Was %d,%d to %d,%d\n", ungappedExtension->start.queryOffset, ungappedExtension->start.subjectOffset,
//    	ungappedExtension->end.queryOffset, ungappedExtension->end.subjectOffset);
//    printf("seed=%d,%d\n", seed.queryOffset, seed.subjectOffset); fflush(stdout);

    // Update the seed to point to byte-packed sequences
	seed.subjectOffset /= 4;

    // Perform dynamic programming for points before the seed
    if (parameters_bytepackStartGap == 0)
    {
    	beforeDpResults = fasterBytepackGappedScoring_dpBeforeSeed(packedSubject, PSSMatrix,
                                                                   seed, dropoff);
    }
    else
    {
    	beforeDpResults = bytepackGappedScoring_dpBeforeSeed(packedSubject, PSSMatrix,
                                                             seed, dropoff);
	}

	// Adjust starting position to point to non-bytepacked subject
	beforeDpResults.best.subjectOffset *= 4;

    // Chop the start off the query and subject so they begin at the seed
    choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
    choppedSubject = packedSubject + seed.subjectOffset;
    choppedSubjectSize = packedSubjectSize - seed.subjectOffset;

    // Perform dynamic programming for points after the seed
    if (parameters_bytepackStartGap == 0)
    {
        afterDpResults = fasterBytepackGappedScoring_dpAfterSeed(choppedSubject, choppedPSSMatrix,
                                                                 dropoff, choppedSubjectSize);
	}
    else
    {
        afterDpResults = bytepackGappedScoring_dpAfterSeed(choppedSubject, choppedPSSMatrix,
                                                           dropoff, choppedSubjectSize);
    }

    // Re-adjust result change due to chopping subject/query and strand adjustment
    afterDpResults.best.queryOffset += seed.queryOffset + strandOffset;
    afterDpResults.best.subjectOffset = (afterDpResults.best.subjectOffset + seed.subjectOffset) * 4;
    beforeDpResults.best.queryOffset += strandOffset;

    // Extend out ends slightly
	beforeDpResults.best.queryOffset -= 2;
	beforeDpResults.best.subjectOffset -= 2;
	afterDpResults.best.queryOffset += 6;
	afterDpResults.best.subjectOffset += 6;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        printf("bytepacked[%d,%d,%d] (seed=%d,%d)\n", beforeDpResults.bestScore, afterDpResults.bestScore,
        PSSMatrix_packedScore[choppedPSSMatrix.bytePackedCodes[0] ^ choppedSubject[0]], seed.queryOffset,
        seed.subjectOffset * 4);
	#endif

    // Determine score by combining score from the two traces, and the match score at
    // the seed position
    score = beforeDpResults.bestScore + afterDpResults.bestScore +
            PSSMatrix_packedScore[choppedPSSMatrix.bytePackedCodes[0] ^ choppedSubject[0]];

	// Use ungapped alignment score and endpoints if better
    if (ungappedExtension->nominalScore > score)
    {
    	score = ungappedExtension->nominalScore;
	}
    else
    {
        // Associate best scoring start and end points with the ungapped extension
        ungappedExtension->start = beforeDpResults.best;
        ungappedExtension->end = afterDpResults.best;
    }

    return score;
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults bytepackGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                    struct coordinate seed, int4 dropoff)
{
    unsigned char *queryPosition, *bestQueryPosition;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, *previousBestRow, *insertQrow, insertS, rowOffset, bestRowToggle = 0;
    int4 match, I, columnCount;
    unsigned char rightOfDropoff;
    struct dpResults dpResults;

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
    if (seed.queryOffset >= bytepackGappedScoring_rowSizes)
    {
        // Free existing rows
        free(bytepackGappedScoring_bestRows[0]);
        free(bytepackGappedScoring_bestRows[1]);
        free(bytepackGappedScoring_insertQrow);

        // Set size to double current needed length
        bytepackGappedScoring_rowSizes = (seed.queryOffset) * 2;

        // Malloc new rows
        bytepackGappedScoring_bestRows[0] = (int4*)global_malloc(sizeof(int4) * bytepackGappedScoring_rowSizes);
        bytepackGappedScoring_bestRows[1] = (int4*)global_malloc(sizeof(int4) * bytepackGappedScoring_rowSizes);
        bytepackGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * bytepackGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject + seed.subjectOffset;
    bestQueryPosition = queryPosition = PSSMatrix.bytePackedCodes + seed.queryOffset;

    // Initialize row pointers
    rowOffset = seed.queryOffset;
   	bestRow = bytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;
    insertQrow = bytepackGappedScoring_insertQrow + rowOffset;

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.bytePackedCodes;
    columnDropoff = PSSMatrix.bytePackedCodes + seed.queryOffset + 1;

    // -----FIRST ROW-----
	// -----SEED POINT----
	*bestRow = 0;
    insertS = -parameters_bytepackOpenGap;
    *insertQrow = -parameters_bytepackOpen4Gap;
    queryPosition--; bestRow--; insertQrow--;

    // -----FIRST CELL-----
    // Only gap open possible
    if (queryPosition >= PSSMatrix.bytePackedCodes)
    {
        *bestRow = insertS;
        insertS = insertS - parameters_bytepackExtendGap;
        *insertQrow = insertS - parameters_bytepackExtend4Gap;
        queryPosition--; bestRow--; insertQrow--;
	}

    // ----- REMAINING CELLS -----
    // For each remaining column in the bottom row, scanning from right-to-left,
    // processing at least 3 more cells
    columnCount = 3;
    while (queryPosition >= PSSMatrix.bytePackedCodes)
    {
        // Only extend gap possible
        *bestRow = insertS;
        insertS = insertS - parameters_bytepackExtendGap;
        *insertQrow = insertS - parameters_bytepackExtend4Gap;

        // If score at current cell is below dropoff
        if (columnCount <= 0 && bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = queryPosition;
            // And stop processing row
            break;
        }

        columnCount--;
        queryPosition--; bestRow--; insertQrow--;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        gappedExtension_printBeforeRow(bytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
                                       rowDropoff, columnDropoff);
	#endif

    // -----REMAINING ROWS-----
    while (subjectPosition > subject && rowDropoff < columnDropoff)
    {
        subjectPosition--;
        queryPosition = columnDropoff - 1;
		rightOfDropoff = 1;

        // Reset row pointers to start of rows and toggle current bestRow
        rowOffset = (queryPosition - PSSMatrix.bytePackedCodes);
        previousBestRow = bytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;
        bestRowToggle = (bestRowToggle + 1) % 2;
        bestRow = bytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;
        insertQrow = bytepackGappedScoring_insertQrow + rowOffset;

        // -----FIRST 4 FAR RIGHT CELLS-----
        columnCount = 4;
        while (columnCount && queryPosition >= rowDropoff)
        {
        	columnCount--;

            // I is the best
            *bestRow = maximum(*insertQrow, insertS);

            // Calculate new I values
            *insertQrow = *bestRow - parameters_bytepackExtend4Gap;
            insertS = *bestRow - parameters_bytepackExtendGap;

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

            queryPosition--; bestRow--; previousBestRow--; insertQrow--;
        }

        // -----CELLS RIGHT OF ROW DROPOFF-----
        while (queryPosition >= rowDropoff)
        {
//        	printf("[[%d]]", previousBestRow[4]);
            // Calculate new M value
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[4];
			I = maximum(insertS, *insertQrow);

            // Match is best
            if (match > I)
            {
            	*bestRow = match;
				if (match - parameters_bytepackStartGap > I)
                {
                	// Match is also best to use for insertion
                	*insertQrow = match - parameters_bytepackOpen4Gap;
                    insertS = match - parameters_bytepackOpenGap;

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
                	// Best to extend insertion
                	*insertQrow = I - parameters_bytepackExtend4Gap;
                    insertS = I - parameters_bytepackExtendGap;
                }
            }
            else
            {
            	// Insertion is best, so also extend it
				*bestRow = I;

                *insertQrow = I - parameters_bytepackExtend4Gap;
                insertS = I - parameters_bytepackExtendGap;
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

            queryPosition--; bestRow--; previousBestRow--; insertQrow--;
        }

        // Record dropoff position
        rowDropoff = queryPosition + 1;

        columnCount = 4;
		while (columnCount && queryPosition >= PSSMatrix.bytePackedCodes)
        {
            // Consider insertion and match
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[4];

            // Match is best
            if (match > insertS)
            {
                *bestRow = match;
                if (match - parameters_bytepackStartGap > insertS)
                {
                    // Match is also best to use for insertion
                    *insertQrow = match - parameters_bytepackOpen4Gap;
                    insertS = match - parameters_bytepackOpenGap;
                }
                else
                {
                    // Best to extend insertion
                    *insertQrow = insertS - parameters_bytepackExtend4Gap;
                    insertS = insertS - parameters_bytepackExtendGap;
                }
            }
            else
            {
                // Insertion is best, so also extend it
                *bestRow = insertS;

                *insertQrow = insertS - parameters_bytepackExtend4Gap;
                insertS = insertS - parameters_bytepackExtendGap;
            }

            if (bestScore <= *bestRow + dropoff)
            {
                // Record dropoff position
                rowDropoff = queryPosition + 1;
			}

            columnCount--;
            queryPosition--; bestRow--; previousBestRow--; insertQrow--;
        }

        // -----CELLS LEFT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow + 1) + dropoff))
        {
            while (queryPosition >= PSSMatrix.bytePackedCodes)
            {
                // Consider insertion only
                *bestRow = insertS;
                *insertQrow = insertS - parameters_bytepackExtend4Gap;
                insertS = insertS - parameters_bytepackExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Stop processing row
                    queryPosition--;
                    break;
                }

                queryPosition--; bestRow--; previousBestRow--; insertQrow--;
            }

            // Record dropoff position
            rowDropoff = queryPosition + 1;
        }

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

	    #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
            gappedExtension_printBeforeRow(bytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
                                           rowDropoff, columnDropoff);
        #endif
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.bytePackedCodes;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}


// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults bytepackGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                   int4 dropoff, int4 subjectLength)
{
    unsigned char *queryPosition, *bestQueryPosition, *queryEnd;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
    int4 bestScore = 0;
    int4 *bestRow, *previousBestRow, *insertQrow, insertS, rowOffset, bestRowToggle = 0;
    int4 match, I, columnCount;
    unsigned char leftOfDropoff;
    int4 queryLength;
    struct dpResults dpResults;

    queryLength = PSSMatrix.length;
    subjectEnd = subject + subjectLength;
    queryEnd = PSSMatrix.bytePackedCodes + queryLength;

//    printf("[queryLength=%d subjectLength=%d]\n", queryLength, subjectLength);

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (queryLength >= bytepackGappedScoring_rowSizes)
    {
        // Free existing rows
        free(bytepackGappedScoring_bestRows[0]);
        free(bytepackGappedScoring_bestRows[1]);
        free(bytepackGappedScoring_insertQrow);

        // Set size to double current needed length
        bytepackGappedScoring_rowSizes = queryLength * 2;

        // Malloc new rows
        bytepackGappedScoring_bestRows[0] = (int4*)global_malloc(sizeof(int4) * bytepackGappedScoring_rowSizes);
        bytepackGappedScoring_bestRows[1] = (int4*)global_malloc(sizeof(int4) * bytepackGappedScoring_rowSizes);
        bytepackGappedScoring_insertQrow = (int4*)global_malloc(sizeof(int4) * bytepackGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject;
    bestQueryPosition = queryPosition = PSSMatrix.bytePackedCodes;

    // Initialize row pointers
   	bestRow = bytepackGappedScoring_bestRows[bestRowToggle];
    insertQrow = bytepackGappedScoring_insertQrow;

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.bytePackedCodes + queryLength - 1;
    columnDropoff = PSSMatrix.bytePackedCodes - 1;

    // -----FIRST ROW-----
	// -----SEED POINT----
	*bestRow = 0;
    insertS = -parameters_bytepackOpenGap;
    *insertQrow = -parameters_bytepackOpen4Gap;
    queryPosition++; bestRow++; insertQrow++;

    // -----FIRST CELL-----
    // Only gap open possible
	if (queryPosition < queryEnd)
    {
        *bestRow = insertS;
        insertS = insertS - parameters_bytepackExtendGap;
        *insertQrow = insertS - parameters_bytepackExtend4Gap;
        queryPosition++; bestRow++; insertQrow++;
	}

    // ----- REMAINING CELLS -----
    // For each remaining column in the bottom row, scanning from right-to-left
    // processing at least 3 more cells
    columnCount = 3;
    while (queryPosition < queryEnd)
    {
        // Only extend gap possible
        *bestRow = insertS;
        insertS = insertS - parameters_bytepackExtendGap;
        *insertQrow = insertS - parameters_bytepackExtend4Gap;

        // If score at current cell is below dropoff
        if (columnCount <= 0 && bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = queryPosition;
            // And stop processing row
            break;
        }

        columnCount--;
        queryPosition++; bestRow++; insertQrow++;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        gappedExtension_printAfterRow(bytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
                                       rowDropoff, columnDropoff);
	#endif

    subjectPosition++;

    // -----REMAINING ROWS-----
    while (subjectPosition < subjectEnd && rowDropoff > columnDropoff)
    {
        queryPosition = columnDropoff + 1;
		leftOfDropoff = 1;

        // Reset row pointers to start of rows and toggle current bestRow
        rowOffset = (queryPosition - PSSMatrix.bytePackedCodes);
        previousBestRow = bytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;
        bestRowToggle = (bestRowToggle + 1) % 2;
        bestRow = bytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;
        insertQrow = bytepackGappedScoring_insertQrow + rowOffset;

        // -----FIRST 4 FAR LEFT CELLS-----
        columnCount = 4;
        while (columnCount && queryPosition <= rowDropoff)
        {
        	columnCount--;

            // I is the best
            *bestRow = maximum(*insertQrow, insertS);

            // Calculate new I values
            *insertQrow = *bestRow - parameters_bytepackExtend4Gap;
            insertS = *bestRow - parameters_bytepackExtendGap;

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

            queryPosition++; bestRow++; previousBestRow++; insertQrow++;
        }

        // -----CELLS LEFT OF ROW DROPOFF-----
        while (queryPosition <= rowDropoff)
        {
//        	printf("[[%d]]", previousBestRow[-4]);
            // Calculate new M value
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[-4];
			I = maximum(insertS, *insertQrow);

            // Match is best
            if (match > I)
            {
            	*bestRow = match;
				if (match - parameters_bytepackStartGap > I)
                {
                	// Match is also best to use for insertion
                	*insertQrow = match - parameters_bytepackOpen4Gap;
                    insertS = match - parameters_bytepackOpenGap;

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
                	// Best to extend insertion
                	*insertQrow = I - parameters_bytepackExtend4Gap;
                    insertS = I - parameters_bytepackExtendGap;
                }
            }
            else
            {
            	// Insertion is best, so also extend it
				*bestRow = I;

                *insertQrow = I - parameters_bytepackExtend4Gap;
                insertS = I - parameters_bytepackExtendGap;
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

            queryPosition++; bestRow++; previousBestRow++; insertQrow++;
        }

        // Record dropoff position
        rowDropoff = queryPosition - 1;

        // Process the first 4 cells right of dropoff
        columnCount = 4;
        while (columnCount && queryPosition < queryEnd)
        {
            // Consider insertion and match
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[-4];

            // Match is best
            if (match > insertS)
            {
                *bestRow = match;
                if (match - parameters_bytepackStartGap > insertS)
                {
                    // Match is also best to use for insertion
                    *insertQrow = match - parameters_bytepackOpen4Gap;
                    insertS = match - parameters_bytepackOpenGap;
                }
                else
                {
                    // Best to extend insertion
                    *insertQrow = insertS - parameters_bytepackExtend4Gap;
                    insertS = insertS - parameters_bytepackExtendGap;
                }
            }
            else
            {
                // Insertion is best, so also extend it
                *bestRow = insertS;

                *insertQrow = insertS - parameters_bytepackExtend4Gap;
                insertS = insertS - parameters_bytepackExtendGap;
            }

            if (bestScore <= *bestRow + dropoff)
            {
                // Record dropoff position
                rowDropoff = queryPosition - 1;
			}

            columnCount--;
            queryPosition++; bestRow++; previousBestRow++; insertQrow++;
		}

        // -----CELLS RIGHT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow - 1) + dropoff))
        {
            while (queryPosition < queryEnd)
            {
                // Consider insertion only
                *bestRow = insertS;
                *insertQrow = insertS - parameters_bytepackExtend4Gap;
                insertS = insertS - parameters_bytepackExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Stop processing row
                    queryPosition++;
                    break;
                }

                queryPosition++; bestRow++; previousBestRow++; insertQrow++;
            }

            // Record dropoff position
            rowDropoff = queryPosition - 1;
        }

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

	    subjectPosition++;

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
            gappedExtension_printAfterRow(bytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
                                           rowDropoff, columnDropoff);
        #endif
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.bytePackedCodes;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

void bytepackGappedScoring_free()
{
    free(bytepackGappedScoring_bestRows[0]);
    free(bytepackGappedScoring_bestRows[1]);
    free(bytepackGappedScoring_insertQrow);
}

