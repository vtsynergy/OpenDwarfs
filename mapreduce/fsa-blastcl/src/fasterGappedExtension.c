// fasterGappedExtension.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for performing a gapped alignment that records traceback information and uses
// the BLAST dropoff technique. This version is faster but less memory efficient than
// gappedExtension; it requires N x M bytes for traceback array

#include "blast.h"

int4* fasterGappedExtension_matchRow = NULL;
int4* fasterGappedExtension_insertQrow = NULL;
int4* fasterGappedExtension_insertSrow = NULL;
int4 fasterGappedExtension_rowSizes = 0;
unsigned char **fasterGappedExtension_traceback = NULL;
int4 fasterGappedExtension_numRows = 0;
int4 fasterGappedExtension_numColumns = 0;

// Prototypes
struct dpResults fasterGappedExtension_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                              struct coordinate seed, int4 dropoff);
struct dpResults fasterGappedExtension_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                             int4 dropoff, int4 subjectLength);
struct trace fasterGappedExtension_traceBeforeSeed(struct dpResults beforeDpResults, struct coordinate seed);
struct trace fasterGappedExtension_traceAfterSeed(struct dpResults beforeDpResults, int4 queryLength);
struct trace fasterGappedExtension_joinTraces(struct trace beforeTrace, struct trace afterTrace);

// Build a gapped extension with a trace and nominal score from the seed point of an ungapped
// extension using dynamic programming
struct gappedExtension* fasterGappedExtension_build(struct ungappedExtension* ungappedExtension,
                                  struct PSSMatrix PSSMatrix, int4 subjectSize,
                                  unsigned char* subject, int4 dropoff)
{
	struct coordinate seed;
	unsigned char *choppedSubject;
	struct dpResults beforeDpResults, afterDpResults;
	struct trace beforeTrace, afterTrace, trace;
	struct PSSMatrix choppedPSSMatrix;
	int4 choppedSubjectSize;
	struct gappedExtension* gappedExtension;
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

//    printf("Seed=%d,%d Length=%d,%d\n", seed.queryOffset, seed.subjectOffset, PSSMatrix.length, subjectSize);

	beforeDpResults = fasterGappedExtension_dpBeforeSeed(subject, PSSMatrix, seed, dropoff);

	// Trace back and create the trace which specifies the first half of the alignment
	beforeTrace = fasterGappedExtension_traceBeforeSeed(beforeDpResults, seed);

	// Chop the start off the query and subject so they begin at the seed
	choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
	choppedSubject = subject + seed.subjectOffset;
	choppedSubjectSize = subjectSize - (seed.subjectOffset);

	// Perform dynamic programming for points after the seed
	afterDpResults = fasterGappedExtension_dpAfterSeed(choppedSubject, choppedPSSMatrix,
	                                              dropoff, choppedSubjectSize);

	// Trace back to get the trace for the seed onwards
	afterTrace = fasterGappedExtension_traceAfterSeed(afterDpResults, choppedPSSMatrix.length);

	// Join afterTrace to the end of beforeTrace
	trace = fasterGappedExtension_joinTraces(beforeTrace, afterTrace);
	free(afterTrace.traceCodes);

    // Adjust coordinates if extension was performed in the second strand
    afterDpResults.best.queryOffset += strandOffset;
    beforeDpResults.best.queryOffset += strandOffset;
	trace.queryStart += strandOffset;

//	printf("Final trace length=%d\n", trace.length);

	// Create gapped extension
	gappedExtension = (struct gappedExtension*)global_malloc(sizeof(struct gappedExtension));
	gappedExtension->trace = trace;
	gappedExtension->next = NULL;

	// Start of afterTrace is end of the gapped extension, but we need to add seed position
	// to get correct offset
	gappedExtension->queryEnd = seed.queryOffset + afterTrace.queryStart + strandOffset;
	gappedExtension->subjectEnd = seed.subjectOffset + afterTrace.subjectStart;

//	if (dloc == 88197331)
//		printf("final[%d,%d,%d](%d)\n", beforeDpResults.bestScore, afterDpResults.bestScore,
//		choppedPSSMatrix.matrix[0][choppedSubject[0]], seed.queryOffset);

    // Determine score by combining score from the two traces, and the match score at
	// the seed position
	gappedExtension->nominalScore = beforeDpResults.bestScore + afterDpResults.bestScore
	                             + choppedPSSMatrix.matrix[0][choppedSubject[0]];

    // Update ungappedExtension start/end
    ungappedExtension->start.queryOffset = trace.queryStart;
    ungappedExtension->end.queryOffset = gappedExtension->queryEnd;
    ungappedExtension->start.subjectOffset = trace.subjectStart;
    ungappedExtension->end.subjectOffset = gappedExtension->subjectEnd;
	ungappedExtension->nominalScore = gappedExtension->nominalScore;

    return gappedExtension;
}

// Given a gapped extension with a nominal score, calculate the normalized score
// and E-Value
void fasterGappedExtension_score(struct gappedExtension* gappedExtension)
{
	gappedExtension->normalizedScore
		= statistics_gappedNominal2normalized(gappedExtension->nominalScore);

	gappedExtension->eValue
		= statistics_gappedCalculateEvalue(gappedExtension->normalizedScore);
}

// Given a gappedExtension and list of ungappedExtensions, prune the latter to
// remove those which overlap/int4ersect with the gappedExtension
void fasterGappedExtension_prune(struct gappedExtension* gappedExtension,
                           struct ungappedExtension* ungappedExtension)
{
	int4 queryPosition, subjectPosition;
	unsigned char* traceCodes;
	int4 count = 0, offset;
	struct trace trace;
	struct ungappedExtension *currentExtension, *previousExtension;

	// Get trace start point and tracecodes
	trace = gappedExtension->trace;
	queryPosition = trace.queryStart;
	subjectPosition = trace.subjectStart;
	traceCodes = trace.traceCodes;

	// Determine start offset
	offset = subjectPosition - queryPosition;

	// For each position in the trace
	while (count < trace.length)
	{
		// Cycle through ungapped extensions and check if they overlap
		// with this part of the gapped extension
		previousExtension = ungappedExtension;
		currentExtension = ungappedExtension->next;
		while (currentExtension != NULL)
		{
			// If on same offset, within start and end of the ungappedExtension
			if (currentExtension->start.subjectOffset -
                currentExtension->start.queryOffset == offset &&
			    currentExtension->start.queryOffset <= queryPosition &&
				currentExtension->end.queryOffset >= queryPosition)
			{
				// Remove current ungappedExtension from the list
				previousExtension->next = currentExtension->next;
				currentExtension = previousExtension->next;
			}
			else
			{
				previousExtension = currentExtension;
				currentExtension = currentExtension->next;
			}
		}

		// A match
		if (traceCodes[count] == 0)
		{
			queryPosition++;
			subjectPosition++;
		}
		// An insertion wrt. query
		else if (traceCodes[count] == 1)
		{
			subjectPosition++;
			offset++;
		}
		// An insertion wrt. subject
		else
		{
			queryPosition++;
			offset--;
		}
		count++;
	}
}

// Join traces before and after the seed together
struct trace fasterGappedExtension_joinTraces(struct trace beforeTrace, struct trace afterTrace)
{
	struct trace joinedTrace;
	unsigned char* traceCodes;
	int4 count;

	// Joined trace will have length equal to sum of lengths of the two traces,
	// and will start where the beforeTrace starts
	joinedTrace.length = beforeTrace.length + afterTrace.length;
	joinedTrace.queryStart = beforeTrace.queryStart;
	joinedTrace.subjectStart = beforeTrace.subjectStart;
	// Add memory to space already allocated by before traceCodes
	traceCodes = (unsigned char*)global_realloc(beforeTrace.traceCodes, sizeof(unsigned char) *
	             joinedTrace.length);

	// Add after trace codes to end in reverse order
	count = 0;
	while (count < afterTrace.length)
	{
		traceCodes[beforeTrace.length + count] = afterTrace.traceCodes[afterTrace.length - count - 1];
		count++;
	}

	joinedTrace.traceCodes = traceCodes;

	return joinedTrace;
}

// Given the results of dynamic programming (a matrix of trace codes and a highest scoring position in
// the matrix) for finding the START of the alignment, performs the simple operation of finding the path
// from the highest scoring point back to the seed
struct trace fasterGappedExtension_traceBeforeSeed(struct dpResults beforeDpResults, struct coordinate seed)
{
	int4 queryPosition, subjectPosition;
	unsigned char** traceback;
	unsigned char traceCode;
	unsigned char state = 0;
	struct trace trace;
	unsigned char* traceCodes;
	uint4 traceCount = 0;

	traceback = beforeDpResults.traceback;
	trace.queryStart = queryPosition = beforeDpResults.best.queryOffset;
	trace.subjectStart = subjectPosition = beforeDpResults.best.subjectOffset;

	// Declare memory for tracecodes; for maximum possible number of codes that could
	// be generated by this trace
	traceCodes = (unsigned char*)global_malloc(sizeof(unsigned char) *
	             (seed.queryOffset - queryPosition + seed.subjectOffset - subjectPosition));

	while (queryPosition < seed.queryOffset && subjectPosition < seed.subjectOffset)
	{
		// Construct the trace
		traceCodes[traceCount] = state;
		traceCount++;

		traceCode = traceback[queryPosition][subjectPosition];

		// If we got to current cell through a MATCH
		if (state == 0)
		{
			// Move to cell before this one
			queryPosition++;
			subjectPosition++;

			// We are only interested in lowest 2 bits of tracecode
			traceCode = traceCode << 6;
			traceCode = traceCode >> 6;

			// Tracecode determines if we matched or inserted here
			state = traceCode;
		}
		// If we got to current cell through an Ix
		else if (state == 1)
		{
			// Move to cell before this one
			subjectPosition++;

			// We are int4erest in bits 3rd and 4th from right
			traceCode = traceCode << 4;
			traceCode = traceCode >> 6;

			// Tracecode determines if we matched or inserted here
			state = traceCode;
		}
		// If we got to current cell through an Iy
		else if (state == 2)
		{
			// Move to cell before this one
			queryPosition++;

			// We are int4erest in bits 5th and 6th from right
			traceCode = traceCode << 2;
			traceCode = traceCode >> 6;

			// Tracecode determines if we matched or inserted here
			state = traceCode;
		}
	}

	// End trace with insertions needed to get us back to the seed
	// (most likely none will be required)
	while (queryPosition < seed.queryOffset)
	{
		traceCodes[traceCount] = 2;
		traceCount++;
		queryPosition++;
	}

	while (subjectPosition < seed.subjectOffset)
	{
		traceCodes[traceCount] = 1;
		traceCount++;
		subjectPosition++;
	}

	trace.traceCodes = traceCodes;
	trace.length = traceCount;

	return trace;
}

// Given the results of dynamic programming (a matrix of trace codes and a highest scoring position in
// the matrix) for finding the END of the alignment, performs the simple operation of finding the path
// from the highest scoring point back to the seed
struct trace fasterGappedExtension_traceAfterSeed(struct dpResults beforeDpResults, int4 queryLength)
{
	int4 queryPosition, subjectPosition;
	unsigned char** traceback;
	unsigned char traceCode;
	unsigned char state = 0;
	struct trace trace;
	unsigned char* traceCodes;
	uint4 traceCount = 0;

	traceback = beforeDpResults.traceback;
	trace.queryStart = 0;
	trace.subjectStart = 0;

	// Start at the end of the alignment
	queryPosition = beforeDpResults.best.queryOffset;
	subjectPosition = beforeDpResults.best.subjectOffset;

	// Declare memory for tracecodes; for maximum possible number of codes that could
	// be generated by this trace
	traceCodes = (unsigned char*)global_malloc(sizeof(unsigned char) * (queryPosition + subjectPosition));

	while (queryPosition > 0 && subjectPosition > 0)
	{
		traceCode = traceback[queryPosition][subjectPosition];
		// If we got to current cell through a MATCH
		if (state == 0)
		{
			// Move to cell before this one
			queryPosition--;
			subjectPosition--;

			// We are only interested in lowest 2 bits of tracecode
			traceCode = traceCode << 6;
			traceCode = traceCode >> 6;

			// Tracecode determines if we matched or inserted here
			state = traceCode;
		}
		// If we got to current cell through an Ix
		else if (state == 1)
		{
			// Move to cell before this one
			subjectPosition--;

			// We are int4erest in bits 3rd and 4th from right
			traceCode = traceCode << 4;
			traceCode = traceCode >> 6;

			// Tracecode determines if we matched or inserted here
			state = traceCode;
		}
		// If we got to current cell through an Iy
		else if (state == 2)
		{
			// Move to cell before this one
			queryPosition--;

			// We are int4erest in bits 5th and 6th from right
			traceCode = traceCode << 2;
			traceCode = traceCode >> 6;

			// Tracecode determines if we matched or inserted here
			state = traceCode;
		}

		// Construct the trace
		traceCodes[traceCount] = state;
		traceCount++;
	}

	// End trace with insertions needed to get us back to the seed
	// (most likely none will be required)
	while (queryPosition > 0)
	{
		traceCodes[traceCount] = 2;
		traceCount++;
		queryPosition--;
	}

	while (subjectPosition > 0)
	{
		traceCodes[traceCount] = 1;
		traceCount++;
		subjectPosition--;
	}

	trace.traceCodes = traceCodes;
	trace.length = traceCount;
	trace.queryStart = beforeDpResults.best.queryOffset;
	trace.subjectStart = beforeDpResults.best.subjectOffset;

	return trace;
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed
struct dpResults fasterGappedExtension_dpBeforeSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                              struct coordinate seed, int4 dropoff)
{
	int2 **queryPosition, **bestQueryPosition;
	int2 *matrixColumn;
	unsigned char *rowDropoff, *columnDropoff;
	unsigned char *subjectPosition, *bestSubjectPosition;
	unsigned char **tracebackRow, *tracebackColumn;
	int4 bestScore = 0, dropoffThreshold;
	int4 *matchRow, *insertQrow, *insertSrow, rowOffset;
	int4 queryDistance, subjectDistance;
	int4 oldMatch, match, previousOldMatch, previousOldInsertS, previousOldInsertQ;
	int4 previousMatch, previousInsertS;
	struct dpResults dpResults;
	unsigned char rightOfDropoff;

	// Declare processing rows for storing match, insert-subject and insert-query values
	// If current malloced rows aren't big enough
	if (seed.subjectOffset >= fasterGappedExtension_rowSizes)
	{
		// Free existing rows
		free(fasterGappedExtension_matchRow);
		free(fasterGappedExtension_insertQrow);
		free(fasterGappedExtension_insertSrow);
		// Set size to double current needed length
		fasterGappedExtension_rowSizes = (seed.subjectOffset) * 2;
		// Malloc new rows
		fasterGappedExtension_matchRow = (int4*)global_malloc(sizeof(int4) * fasterGappedExtension_rowSizes);
		fasterGappedExtension_insertQrow = (int4*)global_malloc(sizeof(int4) * fasterGappedExtension_rowSizes);
		fasterGappedExtension_insertSrow = (int4*)global_malloc(sizeof(int4) * fasterGappedExtension_rowSizes);
	}

	// Determine lowest score before dropoff
	dropoffThreshold = -dropoff;

//    printf("%d,%d (%d,%d)\n", seed.queryOffset, seed.subjectOffset,
//                              fasterGappedExtension_numRows, fasterGappedExtension_numColumns); fflush(stdout);

    // Increase number of columns in traceback array if neccessary
    if (seed.subjectOffset > fasterGappedExtension_numColumns)
    {
    	// For each existing row
        queryDistance = 0;
        while (queryDistance < fasterGappedExtension_numRows)
        {
            // Increase number of columns
            fasterGappedExtension_traceback[queryDistance]
                = (unsigned char*)global_realloc(fasterGappedExtension_traceback[queryDistance],
                                          sizeof(unsigned char) * (seed.subjectOffset));

            queryDistance++;
        }

        // Update number of columns
        fasterGappedExtension_numColumns = seed.subjectOffset;
	}

    // If more rows are required
    if (seed.queryOffset > fasterGappedExtension_numRows)
    {
        // Increase number of row pointers
        fasterGappedExtension_traceback = (unsigned char**)global_realloc(fasterGappedExtension_traceback,
                                    sizeof(unsigned char*) * (seed.queryOffset));

        // Declare new rows
        while (fasterGappedExtension_numRows < seed.queryOffset)
        {
			fasterGappedExtension_traceback[fasterGappedExtension_numRows]
            	= (unsigned char*)global_malloc(sizeof(unsigned char) * (fasterGappedExtension_numColumns));

            fasterGappedExtension_numRows++;
        }
    }

//    printf("%d,%d (%d,%d) AFTER\n", seed.queryOffset, seed.subjectOffset,
//                              fasterGappedExtension_numRows, fasterGappedExtension_numColumns); fflush(stdout);

    bestSubjectPosition = subjectPosition = subject + seed.subjectOffset - 1;
	bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

	// Initialize row pointers
	rowOffset = (subjectPosition - subject);
	matchRow = fasterGappedExtension_matchRow + rowOffset;
	insertQrow = fasterGappedExtension_insertQrow + rowOffset;
	insertSrow = fasterGappedExtension_insertSrow + rowOffset;

	// Set initial row dropoff and column dropoff
	rowDropoff = subject;
	columnDropoff = subject + seed.subjectOffset;

	// Initialize traceback pointers
	tracebackRow = fasterGappedExtension_traceback + (queryPosition - PSSMatrix.matrix);
	tracebackColumn = *tracebackRow + (subjectPosition - subject);

//    printf("[%d,%d]\n", seed.subjectOffset, (subjectPosition - subject)); fflush(stdout);
//    printf("[%d,%d]\n", seed.queryOffset, (queryPosition - PSSMatrix.matrix)); fflush(stdout);

	// -----FIRST ROW-----

	// Using first column of query matrix
	matrixColumn = *queryPosition;

//    printf("[%d]", matrixColumn); fflush(stdout);
//    printf("[%d]", subjectPosition); fflush(stdout);

	// -----FIRST CELL-----
	// Set M value for bottom-right cell
	match = matrixColumn[*subjectPosition];
	*matchRow = match;

	// Set DUMMY Ix and Iy values, which should never be used
	*insertSrow = constants_gappedExtensionDummyValue;
	*insertQrow = constants_gappedExtensionDummyValue;

	// M came from M
	*tracebackColumn = 0;

	// If this is the best-yet scoring cell
	if (match > bestScore)
	{
		// Update best start cell data
		bestScore = *matchRow;
		dropoffThreshold = bestScore - dropoff;
		bestQueryPosition = queryPosition;
		bestSubjectPosition = subjectPosition;
	}

	// Record match and insertS for this about-to-be-previous cell
	previousMatch = match;
	previousInsertS = *insertSrow;

	subjectDistance = 0;
	subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;

	// ----- REMAINING CELLS -----
	// For each remaining column in the bottom row, scanning from right-to-left
	while (subjectPosition >= subject)
	{
		// Set value for M
		match = matrixColumn[*subjectPosition]
		      - parameters_openGap - subjectDistance * parameters_extendGap;
		*matchRow = match;

		// Set value for Ix
		if (previousInsertS - parameters_extendGap > previousMatch - parameters_openGap)
		{
			*insertSrow = previousInsertS - parameters_extendGap;
			// M came from Ix and Ix came from Ix
			*tracebackColumn = 5;
		}
		else
		{
			*insertSrow = previousMatch - parameters_openGap;
			// M came from Ix and Ix came from M
			*tracebackColumn = 1;
		}

		// Set DUMMY Iy value, which should never be used
		*insertQrow = constants_gappedExtensionDummyValue;

		// If this is the best-yet scoring cell
		if (match > bestScore)
		{
			// Update best start cell data
			bestScore = match;
			dropoffThreshold = bestScore - dropoff;
			bestQueryPosition = queryPosition;
			bestSubjectPosition = subjectPosition;
		}

		// If score at current cell is below dropoff
		if (dropoffThreshold > match &&
		    dropoffThreshold > *insertSrow)
		{
			// Record dropoff position
			rowDropoff = subjectPosition;
			// And stop processing row
			break;
		}

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;
		subjectDistance++;
	}

//	print(fasterGappedExtension_matchRow, subject, rowDropoff, columnDropoff);

	queryDistance = 0;

	// -----REMAINING ROWS-----
	while (queryPosition > PSSMatrix.matrix && rowDropoff < columnDropoff)
	{
		queryPosition--; tracebackRow--;
		subjectPosition = columnDropoff - 1;
		tracebackColumn = *tracebackRow + (subjectPosition - subject);

		// Reset row pointers to start of rows
		rowOffset = (subjectPosition - subject);
		matchRow = fasterGappedExtension_matchRow + rowOffset;
		insertQrow = fasterGappedExtension_insertQrow + rowOffset;
		insertSrow = fasterGappedExtension_insertSrow + rowOffset;

		// Using next column of query matrix
		matrixColumn = *queryPosition;

		// -----FAR RIGHT CELL-----
		// Record some old values
		previousOldMatch = *matchRow;
		previousOldInsertQ = *insertQrow;
		previousOldInsertS = *insertSrow;

		// Set Iy value
		if (*insertQrow - parameters_extendGap > *matchRow - parameters_openGap)
		{
			*insertQrow = *insertQrow - parameters_extendGap;
			// Iy is derived from Iy, M is derived from Iy
			*tracebackColumn = 34;
		}
		else
		{
			*insertQrow = *matchRow - parameters_openGap;
			// Iy is derived from M, M is derived from Iy
			*tracebackColumn = 2;
		}

		// Set DUMMY values for M and Iy, which should never be used
		match = *matchRow = constants_gappedExtensionDummyValue;
		*insertSrow = constants_gappedExtensionDummyValue;

		// If score at current cell is below dropoff
		if (dropoffThreshold > *insertQrow)
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

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;

		// -----CELLS RIGHT OF ROW DROPOFF-----
		while (subjectPosition >= rowDropoff)
		{
			// Remember old M value (for cell below this one)
			oldMatch = *matchRow;

			// Calculate new M value
			if (previousOldMatch >= previousOldInsertQ)
			{
				if (previousOldMatch >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldMatch;
					// M is derived from M
					*tracebackColumn = 0;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}
			else
			{
				if (previousOldInsertQ >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertQ;
					// M is derived from Iy
					*tracebackColumn = 2;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}

			*matchRow = match;

			// Record some old values
			previousOldMatch = oldMatch;
			previousOldInsertQ = *insertQrow;
			previousOldInsertS = *insertSrow;

			// Set new Iy value
			if (oldMatch - parameters_openGap >= *insertQrow - parameters_extendGap)
			{
				*insertQrow = oldMatch - parameters_openGap;
				// Iy is derived from M
				// No change to traceback
			}
			else
			{
				*insertQrow = *insertQrow - parameters_extendGap;
				// Iy is derived from Iy
				*tracebackColumn |= 32;
			}
			// Calculate new Ix
			if (previousMatch - parameters_openGap >= previousInsertS - parameters_extendGap)
			{
				*insertSrow = previousMatch - parameters_openGap;
				// Ix is derived from M
				// No change to traceback
			}
			else
			{
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix is derived from Ix
				*tracebackColumn |= 4;
			}

			// If this is the best-yet scoring cell
			if (match > bestScore)
			{
				// Update best start cell data
				bestScore = match;
				dropoffThreshold = bestScore - dropoff;
				bestQueryPosition = queryPosition;
				bestSubjectPosition = subjectPosition;
			}

			// If score at current cell (and cells to its right) are below dropoff
			if (rightOfDropoff)
			{
				if (dropoffThreshold > match &&
					dropoffThreshold > *insertSrow &&
					dropoffThreshold > *insertQrow)
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

			// Record match and insertS for this about-to-be-previous cell
			previousMatch = match;
			previousInsertS = *insertSrow;

			subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;
		}

		// -----CELLS LEFT OF ROW DROPOFF -----
		if (!(dropoffThreshold > previousMatch &&
			dropoffThreshold > previousInsertS &&
			dropoffThreshold > *(insertQrow + 1)))
		{
			while (subjectPosition >= subject)
			{
				// Set value for Ix
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix came from Ix
				*tracebackColumn = 4;

				// Set DUMMY values for M and Ix, which should never be used
				*matchRow = constants_gappedExtensionDummyValue;
				*insertQrow = constants_gappedExtensionDummyValue;

				// If score at current cell is below dropoff
				if (dropoffThreshold > *insertSrow)
				{
					// Stop processing row
					subjectPosition--;
					break;
				}

				// Record match and insertS for this about-to-be-previous cell
				previousInsertS = *insertSrow;

				subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;
				subjectDistance++;
			}
		}

		// Record dropoff position
		rowDropoff = subjectPosition + 1;

//		print(fasterGappedExtension_matchRow, subject, rowDropoff, columnDropoff);

		queryDistance++;
	}

	dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = fasterGappedExtension_traceback;

	return dpResults;
}

// Perform dynamic programming to explore possible END points and alignments that start at
// the given seed
struct dpResults fasterGappedExtension_dpAfterSeed(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                             int4 dropoff, int4 subjectLength)
{
	int2 **queryPosition, **bestQueryPosition, **queryEnd;
	int2 *matrixColumn;
	unsigned char *rowDropoff, *columnDropoff;
	unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
	unsigned char **tracebackRow, *tracebackColumn;
	int4 bestScore = 0, dropoffThreshold;
	int4 *matchRow, *insertQrow, *insertSrow, rowOffset;
	int4 queryDistance, subjectDistance;
	int4 oldMatch, match, previousOldMatch, previousOldInsertS, previousOldInsertQ;
	int4 previousMatch, previousInsertS;
	struct dpResults dpResults;
	unsigned char leftOfDropoff;
	int4 queryLength;

	queryLength = PSSMatrix.length;
	subjectEnd = subject + subjectLength;
	queryEnd = PSSMatrix.matrix + queryLength;

	// Declare processing rows for storing match, insert-subject and insert-query values
	// If current malloced rows aren't big enough
	if (subjectLength >= fasterGappedExtension_rowSizes)
	{
		// Free existing rows
		free(fasterGappedExtension_matchRow);
		free(fasterGappedExtension_insertQrow);
		free(fasterGappedExtension_insertSrow);
		// Set size to double current needed length
		fasterGappedExtension_rowSizes = subjectLength * 2;
		// Malloc new rows
		fasterGappedExtension_matchRow = (int4*)global_malloc(sizeof(int4) * fasterGappedExtension_rowSizes);
		fasterGappedExtension_insertQrow = (int4*)global_malloc(sizeof(int4) * fasterGappedExtension_rowSizes);
		fasterGappedExtension_insertSrow = (int4*)global_malloc(sizeof(int4) * fasterGappedExtension_rowSizes);
	}

	// Determine lowest score before dropoff
	dropoffThreshold = -dropoff;

//    printf("%d,%d (%d,%d) --\n", queryLength, subjectLength,
//                                 fasterGappedExtension_numRows, fasterGappedExtension_numColumns); fflush(stdout);

    // Increase number of columns in traceback array if neccessary
    if (subjectLength > fasterGappedExtension_numColumns)
    {
    	// For each existing row
        queryDistance = 0;
        while (queryDistance < fasterGappedExtension_numRows)
        {
            // Increase number of columns
            fasterGappedExtension_traceback[queryDistance]
                = (unsigned char*)global_realloc(fasterGappedExtension_traceback[queryDistance],
                                          sizeof(unsigned char) * (subjectLength));

            queryDistance++;
        }

        // Update number of columns
        fasterGappedExtension_numColumns = subjectLength;
	}

    // If more rows are required
    if (queryLength > fasterGappedExtension_numRows)
    {
        // Increase number of row pointers
        fasterGappedExtension_traceback = (unsigned char**)global_realloc(fasterGappedExtension_traceback,
                                    sizeof(unsigned char*) * queryLength);

        // Declare new rows
        while (fasterGappedExtension_numRows < queryLength)
        {
			fasterGappedExtension_traceback[fasterGappedExtension_numRows]
            	= (unsigned char*)global_malloc(sizeof(unsigned char) * fasterGappedExtension_numColumns);

            fasterGappedExtension_numRows++;
        }
    }

    bestSubjectPosition = subjectPosition = subject + 1;
	bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

	// Initialize rows
	matchRow = fasterGappedExtension_matchRow + 1;
	insertQrow = fasterGappedExtension_insertQrow + 1;
	insertSrow = fasterGappedExtension_insertSrow + 1;

	// Set initial row dropoff and column dropoff
	rowDropoff = subject + subjectLength - 1;
	columnDropoff = subject;

	// Initialize traceback pointers
	tracebackRow = fasterGappedExtension_traceback + (queryPosition - PSSMatrix.matrix);
	tracebackColumn = *tracebackRow + (subjectPosition - subject);

	// -----FIRST ROW-----

	// Using first column of the query matrix
	matrixColumn = (*queryPosition);

	// -----FIRST CELL-----
	// Set M value for top-left cell
	match = matrixColumn[*subjectPosition];
	*matchRow = match;
	// Set DUMMY Ix and Iy values, which should never be used
	*insertSrow = constants_gappedExtensionDummyValue;
	*insertQrow = constants_gappedExtensionDummyValue;
	// M came from M
	*tracebackColumn = 0;

	// If this is the best-yet scoring cell
	if (match > bestScore)
	{
		// Update best start cell data
		bestScore = match;
		dropoffThreshold = bestScore - dropoff;
		bestQueryPosition = queryPosition;
		bestSubjectPosition = subjectPosition;
	}

	// Record match and insertS for this about-to-be-previous cell
	previousMatch = match;
	previousInsertS = *insertSrow;

	subjectDistance = 0;
	subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;

	// ----- REMAINING CELLS -----
	// For each remaining columns in the top row, scanning from left-to-right
	while (subjectPosition < subjectEnd)
	{
		// Set value for M
		match = matrixColumn[*subjectPosition]
		      - parameters_openGap - subjectDistance * parameters_extendGap;
		*matchRow = match;

		// Set value for Ix
		if (previousInsertS - parameters_extendGap >
			previousMatch - parameters_openGap)
		{
			*insertSrow = previousInsertS - parameters_extendGap;
			// M came from Ix and Ix came from Ix
			*tracebackColumn = 5;
		}
		else
		{
			*insertSrow = previousMatch - parameters_openGap;
			// M came from Ix and Ix came from M
			*tracebackColumn = 1;
		}

		// Set DUMMY Iy value, which should never be used
		*insertQrow = constants_gappedExtensionDummyValue;

		// If this is the best-yet scoring cell
		if (match > bestScore)
		{
			// Update best start cell data
			bestScore = match;
			dropoffThreshold = bestScore - dropoff;
			bestQueryPosition = queryPosition;
			bestSubjectPosition = subjectPosition;
		}

		// If score at current cell is below dropoff
		if (dropoffThreshold > match &&
		    dropoffThreshold > *insertSrow)
		{
			// Record dropoff position
			rowDropoff = subjectPosition;
			// And stop processing row
			break;
		}

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;
		subjectDistance++;
	}

//    if (dloc==88197331)
//    print2(fasterGappedExtension_matchRow, subject, rowDropoff, columnDropoff);

	queryDistance = 0;
	queryPosition++; tracebackRow++;

	// -----REMAINING ROWS-----
	while (queryPosition < queryEnd && rowDropoff > columnDropoff)
	{
		subjectPosition = columnDropoff + 1;
		tracebackColumn = *tracebackRow + (subjectPosition - subject);

		// Reset rows
		rowOffset = (subjectPosition - subject);
		matchRow = fasterGappedExtension_matchRow + rowOffset;
		insertQrow = fasterGappedExtension_insertQrow + rowOffset;
		insertSrow = fasterGappedExtension_insertSrow + rowOffset;

		// Using next column of the query matrix
		matrixColumn = (*queryPosition);

		// -----FAR LEFT CELL-----
		// Record some old values
		previousOldMatch = *matchRow;
		previousOldInsertQ = *insertQrow;
		previousOldInsertS = *insertSrow;

		// Set Iy value
		if (*insertQrow - parameters_extendGap > *matchRow - parameters_openGap)
		{
			*insertQrow = *insertQrow - parameters_extendGap;
			// Iy is derived from Iy, M is derived from Iy
			*tracebackColumn = 34;
		}
		else
		{
			*insertQrow = *matchRow - parameters_openGap;
			// Iy is derived from M, M is derived from Iy
			*tracebackColumn = 2;
		}

		// Set DUMMY values for M and Iy, which should never be used
		match = *matchRow = constants_gappedExtensionDummyValue;
		*insertSrow = constants_gappedExtensionDummyValue;

		// If score at current cell is below dropoff
		if (dropoffThreshold > *insertQrow)
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

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;

		// -----CELLS LEFT OF ROW DROPOFF-----
		while (subjectPosition <= rowDropoff)
		{
			// Remember old M value (for cell below this one)
			oldMatch = *matchRow;

			// Calculate new M value
			if (previousOldMatch >= previousOldInsertQ)
			{
				if (previousOldMatch >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldMatch;
					// M is derived from M
					*tracebackColumn = 0;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}
			else
			{
				if (previousOldInsertQ >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertQ;
					// M is derived from Iy
					*tracebackColumn = 2;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}
			
			*matchRow = match;

			// Record some old values
			previousOldMatch = oldMatch;
			previousOldInsertQ = *insertQrow;
			previousOldInsertS = *insertSrow;

			// Set new Iy value
			if (oldMatch - parameters_openGap >= *insertQrow - parameters_extendGap)
			{
				*insertQrow = oldMatch - parameters_openGap;
				// Iy is derived from M
				// No change to traceback
			}
			else
			{
				*insertQrow = *insertQrow - parameters_extendGap;
				// Iy is derived from Iy
				*tracebackColumn |= 32;
			}
			// Calculate new Ix
			if (previousMatch - parameters_openGap >= previousInsertS - parameters_extendGap)
			{
				*insertSrow = previousMatch - parameters_openGap;
				// Ix is derived from M
				// No change to traceback
			}
			else
			{
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix is derived from Ix
				*tracebackColumn |= 4;
			}

			// If this is the best-yet scoring cell
			if (match > bestScore)
			{
				// Update best start cell data
				bestScore = match;
				dropoffThreshold = bestScore - dropoff;
				bestQueryPosition = queryPosition;
				bestSubjectPosition = subjectPosition;
			}

			// If score at current cell (and cells to its left) are below dropoff
			if (leftOfDropoff)
			{
				if (dropoffThreshold > match &&
					dropoffThreshold > *insertSrow &&
					dropoffThreshold > *insertQrow)
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

			// Record match and insertS for this about-to-be-previous cell
			previousMatch = match;
			previousInsertS = *insertSrow;

			subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;
		}

		// -----CELLS RIGHT OF ROW DROPOFF -----
		if (!(dropoffThreshold > previousMatch &&
			dropoffThreshold > previousInsertS &&
			dropoffThreshold > *(insertQrow - 1)))
		{
			while (subjectPosition < subjectEnd)
			{
				// Set value for Ix
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix came from Ix
				*tracebackColumn = 4;

				// Set DUMMY values for M and Ix, which should never be used
				*matchRow = constants_gappedExtensionDummyValue;
				*insertQrow = constants_gappedExtensionDummyValue;

				// If score at current cell is below dropoff
				if (dropoffThreshold > *insertSrow)
				{
					// And stop processing row
					subjectPosition++;
					break;
				}

				// Record insertS for this about-to-be-previous cell
				previousInsertS = *insertSrow;

				subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;
				subjectDistance++;
			}
		}

		// Record dropoff position
		rowDropoff = subjectPosition - 1;

		queryDistance++;
		queryPosition++; tracebackRow++;

//        if (dloc==88197331)
//		print2(fasterGappedExtension_matchRow, subject, rowDropoff, columnDropoff);
	}

	dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = fasterGappedExtension_traceback;

	return dpResults;
}

void fasterGappedExtension_free()
{
	// Free memory used by traceback array
	while (fasterGappedExtension_numRows > 0)
	{
    	fasterGappedExtension_numRows--;
		free(fasterGappedExtension_traceback[fasterGappedExtension_numRows]);
	}
	free(fasterGappedExtension_traceback);

    free(fasterGappedExtension_matchRow);
    free(fasterGappedExtension_insertQrow);
    free(fasterGappedExtension_insertSrow);
}

// Debugging routine
void fasterGappedExtension_printBeforeRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                    unsigned char* columnDropoff)
{
	unsigned char* subjectPosition = subject;
    int4 *oldRow, *best;

	while (subjectPosition < rowDropoff)
	{
		printf("     ");
		subjectPosition++;
		row++;
	}

    best = oldRow = row;

	while (subjectPosition < columnDropoff)
	{
    	if (*row > *best)
        	best = row;
		row++;
		subjectPosition++;
	}

    row = oldRow; subjectPosition = rowDropoff;

    while (subjectPosition < columnDropoff)
	{
    	if (row == best)
			printf("%4d*", *row);
        else
			printf("%4d ", *row);
		row++;
		subjectPosition++;
	}

	printf("\n");
}

// Debugging routine
void fasterGappedExtension_printAfterRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                   unsigned char* columnDropoff)
{
	unsigned char* subjectPosition = subject;
    int4 *oldRow, *best;

	while (subjectPosition < columnDropoff)
	{
		printf("     ");
		subjectPosition++;
		row++;
	}

    best = oldRow = row;

	while (subjectPosition < rowDropoff)
	{
    	if (*row > *best)
        	best = row;
		row++;
		subjectPosition++;
	}

    row = oldRow; subjectPosition = columnDropoff;

    while (subjectPosition < rowDropoff)
	{
    	if (row == best)
			printf("%4d*", *row);
        else
			printf("%4d ", *row);
		row++;
		subjectPosition++;
	}

	printf("\n");
}
