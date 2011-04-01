// tableGappedScoring.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for performing table-driven gapped alignment

#include "blast.h"

int4* tableGappedScoring_bestRow = NULL;
int4 tableGappedScoring_rowSizes = 0;

#define maximum(a,b) ((a > b) ? a : b)

// Prototypes
struct dpResults tableGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                    struct coordinate seed, int4 dropoff);
struct dpResults tableGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                   int4 dropoff, int4 subjectLength);
void tableGappedScoring_buildTable();

#define maximum(a,b) ((a > b) ? a : b)
#define maximum3(a,b,c) (maximum(maximum(a,b),c))

int4 tableGappedScoring_minChange;
int4 tableGappedScoring_maxChange;

uint2* tableGappedScoring_tablePoint4ers = NULL;
uint2* tableGappedScoring_table = NULL;
uint2* tableGappedScoring_forwardMatchVectors;
uint2* tableGappedScoring_backwardMatchVectors;
char* tableGappedScoring_tableResults;

int4 tableGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                                int4 subjectSize, unsigned char* packedSubject, int4 dropoff)
{
    struct coordinate seed;
    unsigned char *choppedSubject;
    struct PSSMatrix choppedPSSMatrix;
    int4 choppedSubjectSize, score, packedSubjectSize;
    struct dpResults beforeDpResults, afterDpResults;

    // Increment counter number of HSP's semi-gapped
	blast_numSemiGapped++;

    dropoff -= parameters_bytepackDropoffDecrease;

    // Calculate packed subject size
    packedSubjectSize = (subjectSize + 3) / 4;

    seed = ungappedExtension->seed;

    // Adjust to point to start of byte
    seed.queryOffset -= 2;
    seed.subjectOffset -= 2;

    if (tableGappedScoring_table == NULL)
    	tableGappedScoring_buildTable();

    // Update the seed to point to byte-packed sequences
	seed.subjectOffset /= 4;

//    printf("Seed=%d,%d\n", seed.queryOffset, seed.subjectOffset); fflush(stdout);

    // Perform dynamic programming for pointsinitialTable before the seed
    beforeDpResults = tableGappedScoring_dpBeforeSeed(packedSubject, PSSMatrix, seed, dropoff);

//    printf("Score=%d\n", beforeDpResults.bestScore);

	// Adjust starting position to point to non-bytepacked subject
	beforeDpResults.best.subjectOffset *= 4;

    // Chop the start off the query and subject so they begin at the seed
    choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset - 1);
    choppedSubject = packedSubject + seed.subjectOffset;
    choppedSubjectSize = packedSubjectSize - seed.subjectOffset;

//    printf("[[%d]]", PSSMatrix.queryCodes[0]);
//    encoding_printLetters(PSSMatrix.xorCodes[0], 4);

    // Perform dynamic programming for points after the seed
    afterDpResults = tableGappedScoring_dpAfterSeed(choppedSubject, choppedPSSMatrix,
                                                    dropoff, choppedSubjectSize);

    // Re-adjust result change due to chopping subject/query and byte packing
    afterDpResults.best.queryOffset += seed.queryOffset - 1;
    afterDpResults.best.subjectOffset = (afterDpResults.best.subjectOffset + seed.subjectOffset) * 4;

    // Extend out ends slightly
	beforeDpResults.best.queryOffset -= 2;
	beforeDpResults.best.subjectOffset -= 2;
	afterDpResults.best.queryOffset += 6;
	afterDpResults.best.subjectOffset += 6;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
        printf("table[%d,%d](seed=%d,%d)\n", beforeDpResults.bestScore, afterDpResults.bestScore,
                                             seed.queryOffset, seed.subjectOffset * 4);
    #endif

    // Determine score by combining score from the two traces, and the match score at
    // the seed position
    score = beforeDpResults.bestScore + afterDpResults.bestScore;

//	printf("Table score=%d/%d\n", score, ungappedExtension->nominalScore);

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

struct tableBlock
{
	int4 carryIn;
    int4 matchVector;
    uint2* tablePoint4er;
};

// Compare two blocks by examining carryOut values for each columnCarry
int4 tableGappedScoring_compareBlocks(const void* block1,
                                     const void* block2)
{
	const struct tableBlock *b1, *b2;
	int4 columnCarry = tableGappedScoring_minChange;

	b1 = (struct tableBlock*)block1;
	b2 = (struct tableBlock*)block2;

    while (columnCarry <= tableGappedScoring_maxChange)
    {
        if (b1->tablePoint4er[columnCarry] > b2->tablePoint4er[columnCarry])
        {
            return -1;
        }
        else if (b1->tablePoint4er[columnCarry] < b2->tablePoint4er[columnCarry])
        {
            return 1;
        }

        columnCarry++;
	}

    return 0;
}

// Build the table used to perform nucleotide alignment where the subject is compressed
void tableGappedScoring_buildTable()
{
	int4 xorResults, matchVector, matches[4], matchCount;
    int4 carryIn, maxCarry, carryPosition, tempCarryIn;
    int4 range;
	int4 carryIns[4], columnCarry, tempColumnCarry;
    int4 carryOuts[4], matchScore, carryOut, tempCarryOut;
    uint2 *initialTable, *tablePosition, *previousTablePosition;
    struct tableBlock* tableBlocks, *tableBlock;
    int4 best, last, blockCount;

    // Calculate maximum and minimum change between adjacent cells
	tableGappedScoring_minChange = -parameters_bytepackExtendGap;
	tableGappedScoring_maxChange = parameters_matchScore + parameters_bytepackExtendGap;

    range = tableGappedScoring_maxChange - tableGappedScoring_minChange + 1;

    maxCarry = range*range*range*range;

    // Declare memory to hold pointers into the table for every possible column carry
    // and matchVector value
	tableBlocks = (struct tableBlock*)global_malloc(sizeof(struct tableBlock) * maxCarry * 16);

    // Declare memory for the table itself
    initialTable = (uint2*)global_malloc(sizeof(int2) * range * 16 * maxCarry);
	tablePosition = initialTable;

    // For each possible carry in value
    carryIn = 0;
    while (carryIn < maxCarry)
    {
    	// Break carryIn down into 4 values
    	tempCarryIn = carryIn;
    	carryPosition = 0;
        while (carryPosition < 4)
        {
			carryIns[carryPosition] = (tempCarryIn % range) + tableGappedScoring_minChange;
            tempCarryIn /= range;
            carryPosition++;
        }

//        printf("Carry in = %d,%d,%d,%d\n", carryIns[0], carryIns[1], carryIns[2], carryIns[3]);

		// For each of the possibly match vectors between a query character and subject 4-char
        matchVector = 0;
        while (matchVector < 16)
        {
        	matchCount = 0;
            while (matchCount < 4)
            {
            	matches[matchCount] = (matchVector >> matchCount) & 1;
                matchCount++;
            }

        	// Record pointer into the table
            tableBlock = tableBlocks + ((carryIn << 4) | matchVector);
			tableBlock->tablePoint4er = tablePosition - tableGappedScoring_minChange;
			tableBlock->carryIn = carryIn;
			tableBlock->matchVector = matchVector;

			tablePosition += range;

//        	printf("Matches = %d,%d,%d,%d\n", matches[0], matches[1], matches[2], matches[3]);

            // For each column carry
			columnCarry = tableGappedScoring_minChange;
            while (columnCarry <= tableGappedScoring_maxChange)
            {
//				printf("matchVector=%d columnCarry=%d\n", matchVector, columnCarry);

                // Calculate carry out values
				tempColumnCarry = columnCarry;

                // carryOut position 0
                if (matches[0]) matchScore = parameters_matchScore;
                	else matchScore = parameters_mismatchScore;
				carryOuts[0] = maximum3(tempColumnCarry - parameters_bytepackExtendGap,
                	carryIns[0] - parameters_bytepackExtendGap, matchScore) - tempColumnCarry;

//				printf("Max(%d,%d,%d}\n", tempColumnCarry - parameters_bytepackExtendGap,
//                                        carryIns[0] - parameters_bytepackExtendGap,
//                                        matchScore);
//				printf("CarryOuts[0] = %d\n", carryOuts[0]);
                tempColumnCarry += carryOuts[0];

                // carryOut positions 1 - 3
                tempCarryIn = carryIns[0];
                carryPosition = 1;
                while (carryPosition < 4)
                {
                    if (matches[carryPosition]) matchScore = parameters_matchScore;
                    	else matchScore = parameters_mismatchScore;

                    carryOuts[carryPosition] = maximum3(tempColumnCarry - parameters_bytepackExtendGap,
                        tempCarryIn + carryIns[carryPosition] - parameters_bytepackExtendGap,
                        tempCarryIn + matchScore) - tempColumnCarry;

/*				printf("Max(%d,%d,%d} - %d\n", tempColumnCarry - parameters_bytepackExtendGap,
                                               tempCarryIn + carryIns[carryPosition] - parameters_bytepackExtendGap,
                                               tempCarryIn + matchScore,
                                               tempColumnCarry);
					printf("CarryOuts[%d] = %d\n", carryPosition, carryOuts[carryPosition]);*/
                    tempColumnCarry += carryOuts[carryPosition];
                    tempCarryIn += carryIns[carryPosition];

                    carryPosition++;
                }

                // Calculate the final carryOut value from 4 values
				carryOut = 0;
                carryPosition = 4;
                while (carryPosition > 0)
                {
                	carryPosition--;

                    carryOut *= range;
                    carryOut += carryOuts[carryPosition] - tableGappedScoring_minChange;

//                printf("* %d + %d Co=%d\n", range, carryOuts[carryPosition] - tableGappedScoring_minChange, carryOut);
                }

//                if (carryOut > 2000)
//                	printf("CarryOut=%d\n", carryOut);
                tableBlock->tablePoint4er[columnCarry] = carryOut;

                columnCarry++;
            }
        	matchVector++;
        }

        carryIn++;
    }

    // Build table with results for each carryOut
    tableGappedScoring_tableResults = (char*)global_malloc(sizeof(char) * maxCarry);

    // For each possible carry out value
    carryOut = 0;
    while (carryOut < maxCarry)
    {
    	// Break carryIn down into 4 values
    	tempCarryOut = carryOut;
    	carryPosition = 0;
        while (carryPosition < 4)
        {
			carryOuts[carryPosition] = (tempCarryOut % range) + tableGappedScoring_minChange;
            tempCarryOut /= range;
            carryPosition++;
        }

//        printf("Carry out = %d,%d,%d,%d\n", carryOuts[0], carryOuts[1], carryOuts[2], carryOuts[3]);

		// Calculate best and last
        best = 0;
        last = 0;
        carryPosition = 0;
        while (carryPosition < 4)
        {
//        	printf("(%d)", carryOuts[carryPosition]);
			last += carryOuts[carryPosition];
            if (last > best)
            	best = last;
        	carryPosition++;
        }

        // Record best and last
//        printf("best=%d last=%d\n", best, last);
		tableGappedScoring_tableResults[carryOut] = last;

        carryOut++;
	}

    // Sort the table blocks
	qsort(tableBlocks, maxCarry * 16, sizeof(struct tableBlock), tableGappedScoring_compareBlocks);

    // Build the final table
	tableGappedScoring_tablePoint4ers = (uint2*)global_malloc(sizeof(int2) * maxCarry * 16);
    tableGappedScoring_table = (uint2*)global_malloc(sizeof(int2) * range * 16 * maxCarry);
	tablePosition = tableGappedScoring_table;

    // Scan through sorted blocks and look for duplicates
    blockCount = 0;
    previousTablePosition = NULL;
    while (blockCount < maxCarry * 16)
    {
	    tableBlock = tableBlocks + blockCount;
		carryIn = tableBlock->carryIn;
        matchVector = tableBlock->matchVector;

        // If the block is same as the previous one
        if (previousTablePosition &&
            tableGappedScoring_compareBlocks(tableBlock, tableBlocks + blockCount - 1) == 0)
        {
        	// Reuse the block
			tableGappedScoring_tablePoint4ers[(matchVector * maxCarry + carryIn)]
            	= previousTablePosition - tableGappedScoring_table;
        }
        else
        {
        	// Record a new block
            tableGappedScoring_tablePoint4ers[matchVector * maxCarry + carryIn]
                = tablePosition - tableGappedScoring_minChange - tableGappedScoring_table;

            // Copy value into the block
            columnCarry = tableGappedScoring_minChange;
            while (columnCarry <= tableGappedScoring_maxChange)
            {
                *tablePosition = tableBlock->tablePoint4er[columnCarry];
                columnCarry++;
                tablePosition++;
            }
        }

        previousTablePosition = tableGappedScoring_tablePoint4ers[matchVector * maxCarry + carryIn]
                              + tableGappedScoring_table;
    	blockCount++;
    }

//    printf("Size = %d\n", tablePosition - tableGappedScoring_table);

    // Initialize table for converting xor results to 4-bit match/mismatch vector
    tableGappedScoring_backwardMatchVectors = (uint2*)global_malloc(sizeof(int2) * 256);
    tableGappedScoring_forwardMatchVectors = (uint2*)global_malloc(sizeof(int2) * 256);

    // For each possible result of xoring two byte-packed sequences
	xorResults = 0;
    while (xorResults < 256)
    {
    	// Match vectors when performing dynamic programming backwards (before the seed)
		matchVector = 0;
        if (!(xorResults & 0x3)) matchVector |= 1;
        if (!((xorResults >> 2) & 0x3)) matchVector |= 2;
        if (!((xorResults >> 4) & 0x3)) matchVector |= 4;
        if (!((xorResults >> 6) & 0x3)) matchVector |= 8;

		tableGappedScoring_backwardMatchVectors[xorResults] = matchVector * maxCarry;

    	// Match vectors when performing dynamic programming forwards (after the seed)
		matchVector = 0;
        if (!(xorResults & 0x3)) matchVector |= 8;
        if (!((xorResults >> 2) & 0x3)) matchVector |= 4;
        if (!((xorResults >> 4) & 0x3)) matchVector |= 2;
        if (!((xorResults >> 6) & 0x3)) matchVector |= 1;

		tableGappedScoring_forwardMatchVectors[xorResults] = matchVector * maxCarry;


    	xorResults++;
    }
/*
    carryIn = 0;
    while (carryIn < maxCarry)
    {
        matchVector = 0;
        while (matchVector < 16)
        {
            columnCarry = tableGappedScoring_minChange;
            while (columnCarry <= tableGappedScoring_maxChange)
            {
            	if (tableGappedScoring_tablePoint4ers[(carryIn << 4) | matchVector][columnCarry] !=
                    tableGappedScoring_tablePoint4ers2[(carryIn << 4) | matchVector][columnCarry])
                {
                	printf("Error!\n");
				}
                columnCarry++;
            }
            matchVector++;
        }
        carryIn++;
    }
*/

//    tableGappedScoring_testTable();
}

void tableGappedScoring_printCarry(int4 carry)
{
	int4 carryPosition, range;

    range = tableGappedScoring_maxChange - tableGappedScoring_minChange + 1;

    // Break carry down into 4 values
    carryPosition = 0;
    while (carryPosition < 4)
    {
        printf("%d ", ((carry % range) + tableGappedScoring_minChange));
        carry /= range;
        carryPosition++;
    }
    printf("\n");
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults tableGappedScoring_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                 struct coordinate seed, int4 dropoff)
{
    unsigned char *queryPosition, *bestQueryPosition;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, carry, previousBest;
    int4 best;
    unsigned char rightOfDropoff;
    struct dpResults dpResults;
    uint2 matchVector;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (seed.queryOffset >= tableGappedScoring_rowSizes)
    {
        // Free existing rows
        free(tableGappedScoring_bestRow);

        // Set size to double current needed length
        tableGappedScoring_rowSizes = (seed.queryOffset) * 2;

        // Malloc new rows
        tableGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * tableGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject + seed.subjectOffset;
    bestQueryPosition = queryPosition = PSSMatrix.xorCodes + seed.queryOffset;

    // Initialize row pointers
   	bestRow = tableGappedScoring_bestRow + seed.queryOffset;

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.xorCodes;
    columnDropoff = PSSMatrix.xorCodes + seed.queryOffset;

    // -----FIRST ROW-----
	// -----SEED POINT----
	*bestRow = best = 0;
    queryPosition--; bestRow--;

    // ----- REMAINING CELLS -----
    while (queryPosition >= PSSMatrix.xorCodes)
    {
        // Only insertion across possible
        *bestRow = best = best - parameters_bytepackExtendGap;

        // Stop when score drops too low
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position and stop processing row
            rowDropoff = queryPosition;
            break;
        }

        queryPosition--; bestRow--;
	}

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
    gappedExtension_printBeforeRow(tableGappedScoring_bestRow, PSSMatrix.xorCodes,
                                   rowDropoff, columnDropoff);
    #endif

    // -----REMAINING ROWS-----
    while (subjectPosition > subject && rowDropoff < columnDropoff - 1)
    {
        subjectPosition--;
        queryPosition = columnDropoff;

        // Reset row pointer
        bestRow = tableGappedScoring_bestRow + (queryPosition - PSSMatrix.xorCodes);

        // Far right cell
		carry = 0;
        previousBest = *bestRow;
		*bestRow -= parameters_bytepackExtendGap * 4;

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

        queryPosition--; bestRow--;

        // -----CELLS RIGHT OF ROW DROPOFF -----
        while (queryPosition >= rowDropoff)
        {
        	// Get the match vector for query/subect at this position
			matchVector = tableGappedScoring_backwardMatchVectors[*queryPosition ^ *subjectPosition];

/*            printf("[");
            encoding_printLetters(*queryPosition, 4);
            printf(",");
            encoding_printLetters(*subjectPosition, 4);
			printf("]");
            printf("(%d)", matchVector);*/

//            printf("best=[%d,%d]\n", *bestRow, previousBest);
//            printf("columnCarry=%d matchVector=%d carryIn=(%d) ", (*bestRow - previousBest), matchVector, carry);
//			tableGappedScoring_printCarry(carry); fflush(stdout);

//			printf("[%d]", *bestRow - previousBest); fflush(stdout);

			// Lookup value of carry-out in table
            carry = (tableGappedScoring_table + tableGappedScoring_tablePoint4ers[matchVector + carry])
                    [(*bestRow - previousBest)];

//            printf("carryOut=(%d) ", carry);
//			tableGappedScoring_printCarry(carry); fflush(stdout);

			// Use carry-out to update best score for this column
            previousBest = *bestRow;
			*bestRow += tableGappedScoring_tableResults[carry];

            // Check for best scoring cell yet
            if (*bestRow > bestScore)
            {
				bestQueryPosition = queryPosition;
				bestSubjectPosition = subjectPosition;
                bestScore = *bestRow;
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

            queryPosition--; bestRow--;
        }

        // -----CELLS LEFT OF ROW DROPOFF -----
        best = *(bestRow + 1);
        if (!(bestScore > best + dropoff))
        {
            while (queryPosition >= PSSMatrix.xorCodes)
            {
//            	printf("{%d}", best);
                // Only insertion across possible
                *bestRow = best = best - parameters_bytepackExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Stop processing row
                    queryPosition--;
                    break;
                }

                queryPosition--; bestRow--;
			}

            // Record dropoff position
            rowDropoff = queryPosition + 1;
        }

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
            gappedExtension_printBeforeRow(tableGappedScoring_bestRow, PSSMatrix.xorCodes,
                                           rowDropoff, columnDropoff);
        #endif
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.xorCodes;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
struct dpResults tableGappedScoring_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                                int4 dropoff, int4 subjectLength)
{
    unsigned char *queryPosition, *bestQueryPosition, *queryEnd;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
    int4 bestScore = 0;
    int4 *bestRow, carry, previousBest;
    int4 best;
    unsigned char leftOfDropoff;
    int4 queryLength;
    struct dpResults dpResults;
    uint2 matchVector;

    queryLength = PSSMatrix.length;
    subjectEnd = subject + subjectLength;
    queryEnd = PSSMatrix.xorCodes + queryLength;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (queryLength >= tableGappedScoring_rowSizes)
    {
        // Free existing rows
        free(tableGappedScoring_bestRow);

        // Set size to double current needed length
        tableGappedScoring_rowSizes = queryLength * 2;

        // Malloc new rows
        tableGappedScoring_bestRow = (int4*)global_malloc(sizeof(int4) * tableGappedScoring_rowSizes);
    }

    bestSubjectPosition = subjectPosition = subject;
    bestQueryPosition = queryPosition = PSSMatrix.xorCodes;

    // Initialize row pointers
   	bestRow = tableGappedScoring_bestRow;

    // Set initial row dropoff and column dropoff
    rowDropoff = PSSMatrix.xorCodes + queryLength;
    columnDropoff = PSSMatrix.xorCodes;

    // -----FIRST ROW-----
	// -----SEED POINT----
	*bestRow = best = 0;
    queryPosition++; bestRow++;

    // ----- REMAINING CELLS -----
    while (queryPosition < queryEnd)
    {
        // Only insertion across possible
        *bestRow = best = best - parameters_bytepackExtendGap;

        // Stop when score drops too low
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position and stop processing row
            rowDropoff = queryPosition;
            break;
        }

        queryPosition++; bestRow++;
	}

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
    	gappedExtension_printAfterRow(tableGappedScoring_bestRow, PSSMatrix.xorCodes,
	                                  rowDropoff, columnDropoff);
    #endif

//	subjectPosition++;

    // -----REMAINING ROWS-----
    while (subjectPosition < subjectEnd && rowDropoff > columnDropoff + 1)
    {
        queryPosition = columnDropoff;

        // Reset row pointer
        bestRow = tableGappedScoring_bestRow + (queryPosition - PSSMatrix.xorCodes);

        // First cell
		carry = 0;
        previousBest = *bestRow;
        *bestRow -= parameters_bytepackExtendGap * 4;

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

        queryPosition++; bestRow++;

        // -----CELLS BEFORE ROW DROPOFF -----
        while (queryPosition < rowDropoff)
        {
        	// Get the match vector for query/subect at this position
			matchVector = tableGappedScoring_forwardMatchVectors[*queryPosition ^ *subjectPosition];

/*            printf("[");
            encoding_printLetters(*queryPosition, 4);
            printf(",");
            encoding_printLetters(*subjectPosition, 4);
			printf("]");
            printf("(%d)", matchVector);*/

//            printf("best=[%d,%d]\n", *bestRow, previousBest);
//            printf("columnCarry=%d matchVector=%d carryIn=(%d) ", (*bestRow - previousBest), matchVector, carry);
//			tableGappedScoring_printCarry(carry);

			// Lookup value of carry-out in table
            carry = (tableGappedScoring_table + tableGappedScoring_tablePoint4ers[matchVector + carry])
                    [(*bestRow - previousBest)];

//            printf("last=%d + %d carryOut=(%d) ", tableGappedScoring_tableResults[carry], *bestRow, carry);
//			tableGappedScoring_printCarry(carry); fflush(stdout);

			// Use carry-out to update best score for this column
            previousBest = *bestRow;
			*bestRow += tableGappedScoring_tableResults[carry];

            // Check for best scoring cell yet
            if (*bestRow > bestScore)
            {
				bestQueryPosition = queryPosition;
				bestSubjectPosition = subjectPosition;
                bestScore = *bestRow;
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

            queryPosition++; bestRow++;
        }

        // -----CELLS AFTER ROW DROPOFF -----
        best = *(bestRow - 1);
        if (!(bestScore > best + dropoff))
        {
            while (queryPosition < queryEnd)
            {
//            	printf("{%d}", best);
                // Only insertion across possible
                *bestRow = best = best - parameters_bytepackExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Stop processing row
                    queryPosition++;
                    break;
                }

                queryPosition++; bestRow++;
			}
            // Record dropoff position
            rowDropoff = queryPosition - 1;
        }

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
            gappedExtension_printAfterRow(tableGappedScoring_bestRow, PSSMatrix.xorCodes,
                                          rowDropoff, columnDropoff);
        #endif

        subjectPosition++;
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.xorCodes;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}


void tableGappedScoring_free()
{
}


