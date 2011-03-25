// fasterBytepackGappedScoring.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform gapped alignment using byte-packed alignment technique
// with non-affine gap costs

#include "blast.h"

int4* fasterBytepackGappedScoring_bestRows[2] = {NULL,NULL};
int4 fasterBytepackGappedScoring_rowSizes = 0;

#define maximum(a,b) ((a > b) ? a : b)

// Prototypes
struct dpResults fasterBytepackGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                    struct coordinate seed, int4 dropoff);
struct dpResults fasterBytepackGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                   int4 dropoff, int4 subjectLength);

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults fasterBytepackGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                    struct coordinate seed, int4 dropoff)
{
    unsigned char *queryPosition, *bestQueryPosition;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, *previousBestRow, insertS, rowOffset, bestRowToggle = 0;
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
    if (seed.queryOffset >= fasterBytepackGappedScoring_rowSizes)
    {
        // Free existing rows
        free(fasterBytepackGappedScoring_bestRows[0]);
        free(fasterBytepackGappedScoring_bestRows[1]);

        // Set size to double current needed length
        fasterBytepackGappedScoring_rowSizes = (seed.queryOffset) * 2;

        // Malloc new rows
        fasterBytepackGappedScoring_bestRows[0] = (int4*)global_malloc(sizeof(int4) * fasterBytepackGappedScoring_rowSizes);
        fasterBytepackGappedScoring_bestRows[1] = (int4*)global_malloc(sizeof(int4) * fasterBytepackGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject + seed.subjectOffset;
    bestQueryPosition = queryPosition = PSSMatrix.bytePackedCodes + seed.queryOffset;

    // Initialize row pointers
    rowOffset = seed.queryOffset;
   	bestRow = fasterBytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.bytePackedCodes;
    columnDropoff = PSSMatrix.bytePackedCodes + seed.queryOffset + 1;

    // -----FIRST ROW-----
	// -----SEED POINT----
	insertS = *bestRow = 0;
    queryPosition--; bestRow--;

    // -----FIRST CELL-----
    // Only gap open possible
    if (queryPosition >= PSSMatrix.bytePackedCodes)
    {
        *bestRow = insertS = insertS - parameters_bytepackExtendGap;
        queryPosition--; bestRow--;
	}

    // ----- REMAINING CELLS -----
    // For each remaining column in the bottom row, scanning from right-to-left,
    // processing at least 3 more cells
    columnCount = 3;
    while (queryPosition >= PSSMatrix.bytePackedCodes)
    {
        // Only extend gap possible
        *bestRow = insertS = insertS - parameters_extendGap;

        // If score at current cell is below dropoff
        if (columnCount <= 0 && bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = queryPosition;
            // And stop processing row
            break;
        }

        columnCount--;
        queryPosition--; bestRow--;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        gappedExtension_printBeforeRow(fasterBytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
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
        previousBestRow = fasterBytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;
        bestRowToggle = (bestRowToggle + 1) % 2;
        bestRow = fasterBytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;

        // -----FIRST 4 FAR RIGHT CELLS-----
        columnCount = 4;
        while (columnCount && queryPosition >= rowDropoff)
        {
        	columnCount--;

            // I is the best
            *bestRow = maximum(*previousBestRow - parameters_bytepackExtend4Gap,
                               insertS - parameters_bytepackExtendGap);
			insertS = *bestRow;

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

            queryPosition--; bestRow--; previousBestRow--;
        }

        // -----CELLS RIGHT OF ROW DROPOFF-----
        while (queryPosition >= rowDropoff)
        {
//        	printf("[[%d]]", previousBestRow[4]);
            // Calculate new M value
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[4];

            I = maximum(*previousBestRow - parameters_bytepackExtend4Gap,
                        insertS - parameters_bytepackExtendGap);

            // Match is best
            if (match > I)
            {
            	insertS = *bestRow = match;

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
            	// Insertion is best
				insertS = *bestRow = I;
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

            queryPosition--; bestRow--; previousBestRow--;
        }

        // Record dropoff position
        rowDropoff = queryPosition + 1;

        columnCount = 4;
		while (columnCount && queryPosition >= PSSMatrix.bytePackedCodes)
        {
            // Consider insertion and match
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[4];
            I = insertS - parameters_bytepackExtendGap;

            // Match is best
            if (match > I)
            {
                insertS = *bestRow = match;

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
            	// Insertion is best
				insertS = *bestRow = I;
            }

            if (bestScore <= *bestRow + dropoff)
            {
                // Record dropoff position
                rowDropoff = queryPosition + 1;
			}

            columnCount--;
            queryPosition--; bestRow--; previousBestRow--;
        }

        // -----CELLS LEFT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow + 1) + dropoff))
        {
            while (queryPosition >= PSSMatrix.bytePackedCodes)
            {
                // Consider insertion only
                *bestRow = insertS = insertS - parameters_bytepackExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Stop processing row
                    queryPosition--;
                    break;
                }

                queryPosition--; bestRow--; previousBestRow--;
            }

            // Record dropoff position
            rowDropoff = queryPosition + 1;
        }

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

	    #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
            gappedExtension_printBeforeRow(fasterBytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
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
struct dpResults fasterBytepackGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                   int4 dropoff, int4 subjectLength)
{
    unsigned char *queryPosition, *bestQueryPosition, *queryEnd;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
    int4 bestScore = 0;
    int4 *bestRow, *previousBestRow, insertS, rowOffset, bestRowToggle = 0;
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
    if (queryLength >= fasterBytepackGappedScoring_rowSizes)
    {
        // Free existing rows
        free(fasterBytepackGappedScoring_bestRows[0]);
        free(fasterBytepackGappedScoring_bestRows[1]);

        // Set size to double current needed length
        fasterBytepackGappedScoring_rowSizes = queryLength * 2;

        // Malloc new rows
        fasterBytepackGappedScoring_bestRows[0] = (int4*)global_malloc(sizeof(int4) * fasterBytepackGappedScoring_rowSizes);
        fasterBytepackGappedScoring_bestRows[1] = (int4*)global_malloc(sizeof(int4) * fasterBytepackGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject;
    bestQueryPosition = queryPosition = PSSMatrix.bytePackedCodes;

    // Initialize row pointers
   	bestRow = fasterBytepackGappedScoring_bestRows[bestRowToggle];

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.bytePackedCodes + queryLength - 1;
    columnDropoff = PSSMatrix.bytePackedCodes - 1;

    // -----FIRST ROW-----
	// -----SEED POINT----
	insertS = *bestRow = 0;
    queryPosition++; bestRow++;

    // -----FIRST CELL-----
    // Only gap open possible
	if (queryPosition < queryEnd)
    {
        *bestRow = insertS = insertS - parameters_bytepackExtendGap;
        queryPosition++; bestRow++;
	}

    // ----- REMAINING CELLS -----
    // For each remaining column in the bottom row, scanning from right-to-left
    // processing at least 3 more cells
    columnCount = 3;
    while (queryPosition < queryEnd)
    {
        // Only extend gap possible
        *bestRow = insertS = insertS - parameters_bytepackExtendGap;

        // If score at current cell is below dropoff
        if (columnCount <= 0 && bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = queryPosition;
            // And stop processing row
            break;
        }

        columnCount--;
        queryPosition++; bestRow++;
    }

    // Clear insertS for next row
    insertS = constants_gappedExtensionDummyValue;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        gappedExtension_printAfterRow(fasterBytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
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
        previousBestRow = fasterBytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;
        bestRowToggle = (bestRowToggle + 1) % 2;
        bestRow = fasterBytepackGappedScoring_bestRows[bestRowToggle] + rowOffset;

        // -----FIRST 4 FAR LEFT CELLS-----
        columnCount = 4;
        while (columnCount && queryPosition <= rowDropoff)
        {
        	columnCount--;

            // I is the best
            *bestRow = maximum(*previousBestRow - parameters_bytepackExtend4Gap,
                               insertS - parameters_bytepackExtendGap);
			insertS = *bestRow;

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

            queryPosition++; bestRow++; previousBestRow++;
        }

        // -----CELLS LEFT OF ROW DROPOFF-----
        while (queryPosition <= rowDropoff)
        {
//        	printf("[[%d]]", previousBestRow[-4]);
            // Calculate new M value
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[-4];

            I = maximum(*previousBestRow - parameters_bytepackExtend4Gap,
                        insertS - parameters_bytepackExtendGap);

            // Match is best
            if (match > I)
            {
            	insertS = *bestRow = match;

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
            	// Insertion is best
				insertS = *bestRow = I;
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

            queryPosition++; bestRow++; previousBestRow++;
        }

        // Record dropoff position
        rowDropoff = queryPosition - 1;

        // Process the first 4 cells right of dropoff
        columnCount = 4;
        while (columnCount && queryPosition < queryEnd)
        {
            // Consider insertion and match
            match = PSSMatrix_packedScore[*queryPosition ^ *subjectPosition] + previousBestRow[-4];

            I = insertS - parameters_bytepackExtendGap;

            // Match is best
            if (match > I)
            {
                insertS = *bestRow = match;

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
            	// Insertion is best
				insertS = *bestRow = I;
            }

            if (bestScore <= *bestRow + dropoff)
            {
                // Record dropoff position
                rowDropoff = queryPosition - 1;
			}

            columnCount--;
            queryPosition++; bestRow++; previousBestRow++;
		}

        // -----CELLS RIGHT OF ROW DROPOFF -----
        if (!(bestScore > *(bestRow - 1) + dropoff))
        {
            while (queryPosition < queryEnd)
            {
                // Consider insertion only
                *bestRow = insertS = insertS - parameters_bytepackExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Stop processing row
                    queryPosition++;
                    break;
                }

                queryPosition++; bestRow++; previousBestRow++;
            }

            // Record dropoff position
            rowDropoff = queryPosition - 1;
        }

        // Clear insertS for next row
        insertS = constants_gappedExtensionDummyValue;

	    subjectPosition++;

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
        {
            gappedExtension_printAfterRow(fasterBytepackGappedScoring_bestRows[bestRowToggle], PSSMatrix.bytePackedCodes,
                                           rowDropoff, columnDropoff);
		}
        #endif
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.bytePackedCodes;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

void fasterBytepackGappedScoring_free()
{
    free(fasterBytepackGappedScoring_bestRows[0]);
    free(fasterBytepackGappedScoring_bestRows[1]);
}

