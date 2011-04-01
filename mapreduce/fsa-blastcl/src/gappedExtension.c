// gappedExtension.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for performing a gapped alignment that records traceback information and uses
// the BLAST dropoff technique. This code is slower but more memory efficient than
// fasterGappedExtension, and allows for cases where only part of the subject has been unpacked

#include "blast.h"

int4* gappedExtension_matchRow = NULL;
int4* gappedExtension_insertQrow = NULL;
int4* gappedExtension_insertSrow = NULL;
int4 gappedExtension_rowSizes = 0;

unsigned char **gappedExtension_traceback = NULL;
unsigned char *gappedExtension_tracebackData = NULL;
int4 gappedExtension_numRows = 0;
int4 gappedExtension_tracebackAlloc = 0;

// Prototypes
struct dpResults gappedExtension_dpBeforeSeed(struct PSSMatrix PSSMatrix, int4 dropoff,
       struct coordinate seed, struct unpackRegion* unpackRegion);
struct dpResults gappedExtension_dpAfterSeed(struct PSSMatrix PSSMatrix, int4 dropoff,
       struct unpackRegion* unpackRegion, int4 subjectLength,
       int4 seedSubjectOffset);
struct trace gappedExtension_traceBeforeSeed(struct dpResults beforeDpResults, struct coordinate seed);
struct trace gappedExtension_traceAfterSeed(struct dpResults beforeDpResults, int4 queryLength);
struct trace gappedExtension_joinTraces(struct trace beforeTrace, struct trace afterTrace);

// Build a gapped extension with a trace and nominal score from the seed point of an ungapped
// extension using dynamic programming
struct gappedExtension* gappedExtension_build(struct ungappedExtension* ungappedExtension,
                        struct PSSMatrix PSSMatrix, int4 subjectSize, unsigned char* subject,
                        struct unpackRegion* unpackRegion, int4 dropoff)
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

    beforeDpResults = gappedExtension_dpBeforeSeed(PSSMatrix, dropoff, seed, unpackRegion);

	// Trace back and create the trace which specifies the first half of the alignment
	beforeTrace = gappedExtension_traceBeforeSeed(beforeDpResults, seed);

	// Chop the start off the query and subject so they begin at the seed
	choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
	choppedSubject = subject + seed.subjectOffset;
	choppedSubjectSize = subjectSize - (seed.subjectOffset);

	// Perform dynamic programming for points after the seed
	afterDpResults = gappedExtension_dpAfterSeed(choppedPSSMatrix, dropoff, unpackRegion,
                                                 choppedSubjectSize, seed.subjectOffset);

	// Trace back to get the trace for the seed onwards
	afterTrace = gappedExtension_traceAfterSeed(afterDpResults, choppedPSSMatrix.length);

	// Join afterTrace to the end of beforeTrace
	trace = gappedExtension_joinTraces(beforeTrace, afterTrace);
	free(afterTrace.traceCodes);

    // Adjust coordinates if extension was performed in the second strand
    afterDpResults.best.queryOffset += strandOffset;
    beforeDpResults.best.queryOffset += strandOffset;
	trace.queryStart += strandOffset;

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
//		choppedPSSMatrix.matrix[0][unpackRegion->unpackedSubject[seed.subjectOffset]], seed.queryOffset);

	// Determine score by combining score from the two traces, and the match score at
	// the seed position
	gappedExtension->nominalScore = beforeDpResults.bestScore + afterDpResults.bestScore
	               + choppedPSSMatrix.matrix[0][unpackRegion->unpackedSubject[seed.subjectOffset]];

    // Update ungappedExtension start/end
    ungappedExtension->start.queryOffset = trace.queryStart;
    ungappedExtension->end.queryOffset = gappedExtension->queryEnd;
    ungappedExtension->start.subjectOffset = trace.subjectStart;
    ungappedExtension->end.subjectOffset = gappedExtension->subjectEnd;
	ungappedExtension->nominalScore = gappedExtension->nominalScore;

    #ifdef VERBOSE
    if (parameters_verboseDloc == blast_dloc)
    {
        printf("Gapped Extension from %d,%d to %d,%d score %d\n", trace.queryStart, trace.subjectStart,
               gappedExtension->queryEnd, gappedExtension->subjectEnd, gappedExtension->nominalScore);
    }
    #endif

    return gappedExtension;
}

// Given a gapped extension with a nominal score, calculate the normalized score
// and E-Value
void gappedExtension_score(struct gappedExtension* gappedExtension)
{
	gappedExtension->normalizedScore
		= statistics_gappedNominal2normalized(gappedExtension->nominalScore);

	gappedExtension->eValue
		= statistics_gappedCalculateEvalue(gappedExtension->normalizedScore);
}

// Given a gappedExtension and list of ungappedExtensions, prune the latter to
// remove those which overlap/int4ersect with the gappedExtension
void gappedExtension_prune(struct gappedExtension* gappedExtension,
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
struct trace gappedExtension_joinTraces(struct trace beforeTrace, struct trace afterTrace)
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
struct trace gappedExtension_traceBeforeSeed(struct dpResults beforeDpResults, struct coordinate seed)
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

//        printf("(%p)", traceback[queryPosition]);
//        printf("(%d,%d)", queryPosition, subjectPosition); fflush(stdout);
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
struct trace gappedExtension_traceAfterSeed(struct dpResults beforeDpResults, int4 queryLength)
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
struct dpResults gappedExtension_dpBeforeSeed(struct PSSMatrix PSSMatrix, int4 dropoff,
       struct coordinate seed, struct unpackRegion* unpackRegion)
{
	int2 **queryPosition, **bestQueryPosition;
	int2 *matrixColumn;
	unsigned char *subject, *rowDropoff, *columnDropoff, *minRowDropoff;
	unsigned char *subjectPosition, *bestSubjectPosition;
	unsigned char **tracebackRow, *tracebackColumn;
	int4 bestScore = 0, dropoffThreshold;
	int4 *matchRow, *insertQrow, *insertSrow, rowOffset;
	int4 queryDistance, subjectDistance;
	int4 oldMatch, match, previousOldMatch, previousOldInsertS, previousOldInsertQ, previousInsertQ = 0;
	int4 previousMatch, previousInsertS;
	struct dpResults dpResults;
	unsigned char rightOfDropoff;
    unsigned char *tracebackDataPtr, *newTracebackData;

    subject = unpackRegion->unpackedSubject;

	if (gappedExtension_matchRow == NULL)
    {
        // Declare processing rows for storing match, insert-subject and insert-query values
		gappedExtension_rowSizes = 2 * statistics_gappedFinalNominalDropoff / parameters_extendGap;

		// Malloc new rows
		gappedExtension_matchRow = (int4*)global_malloc(sizeof(int4) * gappedExtension_rowSizes);
		gappedExtension_insertQrow = (int4*)global_malloc(sizeof(int4) * gappedExtension_rowSizes);
		gappedExtension_insertSrow = (int4*)global_malloc(sizeof(int4) * gappedExtension_rowSizes);
    }

	// Determine lowest score before dropoff
	dropoffThreshold = -dropoff;

    // If more rows are required
    if (seed.queryOffset > gappedExtension_numRows)
    {
        // Increase number of row pointers
    	gappedExtension_numRows = seed.queryOffset;
        gappedExtension_traceback = (unsigned char**)global_realloc(gappedExtension_traceback,
                                    sizeof(unsigned char*) * (gappedExtension_numRows));
    }

    // If first time, initialize traceback rows
    if (gappedExtension_tracebackData == NULL)
    {
    	gappedExtension_tracebackAlloc = constants_initialTracebackAlloc;
		gappedExtension_tracebackData = global_malloc(sizeof(char) * gappedExtension_tracebackAlloc);
    }

    tracebackDataPtr = gappedExtension_tracebackData + gappedExtension_tracebackAlloc - 1;

    bestSubjectPosition = subjectPosition = subject + seed.subjectOffset - 1;
	bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

	// Initialize row pointers
	rowOffset = (subjectPosition - subject);
	matchRow = gappedExtension_matchRow + (rowOffset % gappedExtension_rowSizes);
	insertQrow = gappedExtension_insertQrow + (rowOffset % gappedExtension_rowSizes);
	insertSrow = gappedExtension_insertSrow + (rowOffset % gappedExtension_rowSizes);

	// Set initial row dropoff and column dropoff
	rowDropoff = subject;
	columnDropoff = subject + seed.subjectOffset;

    // Calculate minimum row dropoff after this row is processed
    minRowDropoff = columnDropoff - ((dropoff / parameters_extendGap) + 2);
    if (minRowDropoff < subject)
        minRowDropoff = subject;

    // Unpack more of the subject if required
    unpack_extendRegionStart(minRowDropoff - subject, unpackRegion);

    // If the unpacked subject has moved in memory
    if (subject != unpackRegion->unpackedSubject)
    {
        // Update all pointers to point at the new subject
        subjectPosition = (subjectPosition - subject) + unpackRegion->unpackedSubject;
        rowDropoff = (rowDropoff - subject) + unpackRegion->unpackedSubject;
        columnDropoff = (columnDropoff - subject) + unpackRegion->unpackedSubject;
        bestSubjectPosition = (bestSubjectPosition - subject) + unpackRegion->unpackedSubject;
        subject = unpackRegion->unpackedSubject;
    }

    // Initialize traceback pointer for this row
	tracebackRow = gappedExtension_traceback + (queryPosition - PSSMatrix.matrix);
	*tracebackRow = tracebackDataPtr - (subjectPosition - subject);
    tracebackColumn = tracebackDataPtr;
    tracebackDataPtr -= (subjectPosition - minRowDropoff + 1);

	// -----FIRST ROW-----

	// Using first column of query matrix
	matrixColumn = *queryPosition;

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

    // Check for scoring row wrap-around
    if (matchRow < gappedExtension_matchRow)
    {
		matchRow += gappedExtension_rowSizes;
		insertSrow += gappedExtension_rowSizes;
		insertQrow += gappedExtension_rowSizes;
    }

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

        // Check for scoring row wrap-around
        if (matchRow < gappedExtension_matchRow)
        {
            matchRow += gappedExtension_rowSizes;
            insertSrow += gappedExtension_rowSizes;
            insertQrow += gappedExtension_rowSizes;
        }

		subjectDistance++;
	}

//	print(gappedExtension_matchRow, subject, rowDropoff, columnDropoff);

	queryDistance = 0;

	// -----REMAINING ROWS-----
	while (queryPosition > PSSMatrix.matrix && rowDropoff < columnDropoff)
	{
		queryPosition--; tracebackRow--;
		subjectPosition = columnDropoff - 1;

        // Calculate minimum row dropoff after this row is processed
        minRowDropoff = rowDropoff - ((PSSMatrix.highestValue / parameters_extendGap) + 2);
        if (minRowDropoff < subject)
        	minRowDropoff = subject;

        // If not enough space in traceback data for this row, realloc
        if (subjectPosition - minRowDropoff >= tracebackDataPtr - gappedExtension_tracebackData)
        {
			newTracebackData = (unsigned char*)global_malloc(sizeof(char) * gappedExtension_tracebackAlloc * 2);

            // Move existing data to end of new block
            memcpy(newTracebackData + gappedExtension_tracebackAlloc,
                   gappedExtension_tracebackData, gappedExtension_tracebackAlloc);

            // Update traceback data pointer
			tracebackDataPtr = newTracebackData + gappedExtension_tracebackAlloc
                              + (tracebackDataPtr - gappedExtension_tracebackData);

            // Update existing rows to point to new data block
			while (tracebackRow < gappedExtension_traceback + seed.queryOffset - 1)
            {
				tracebackRow++;
            	*tracebackRow = newTracebackData + gappedExtension_tracebackAlloc
                              + (*tracebackRow - gappedExtension_tracebackData);
            }

            // Data block is now double the size
            gappedExtension_tracebackAlloc *= 2;

            free(gappedExtension_tracebackData);
            gappedExtension_tracebackData = newTracebackData;
        }

        // Initialize traceback pointer for this row
        tracebackRow = gappedExtension_traceback + (queryPosition - PSSMatrix.matrix);
        *tracebackRow = tracebackDataPtr - (subjectPosition - subject);
        tracebackColumn = tracebackDataPtr;
        tracebackDataPtr -= (subjectPosition - minRowDropoff + 1);

		// Unpack more of the subject if required
        unpack_extendRegionStart(minRowDropoff - subject, unpackRegion);

        // If the unpacked subject has moved in memory
        if (subject != unpackRegion->unpackedSubject)
        {
        	// Update all pointers to point at the new subject
        	subjectPosition = (subjectPosition - subject) + unpackRegion->unpackedSubject;
			rowDropoff = (rowDropoff - subject) + unpackRegion->unpackedSubject;
			columnDropoff = (columnDropoff - subject) + unpackRegion->unpackedSubject;
            bestSubjectPosition = (bestSubjectPosition - subject) + unpackRegion->unpackedSubject;
            subject = unpackRegion->unpackedSubject;
        }

        // Reset row pointers to start of rows
		rowOffset = (subjectPosition - subject);
        matchRow = gappedExtension_matchRow + (rowOffset % gappedExtension_rowSizes);
        insertQrow = gappedExtension_insertQrow + (rowOffset % gappedExtension_rowSizes);
        insertSrow = gappedExtension_insertSrow + (rowOffset % gappedExtension_rowSizes);

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

        // Check for scoring row wrap-around
        if (matchRow < gappedExtension_matchRow)
        {
            matchRow += gappedExtension_rowSizes;
            insertSrow += gappedExtension_rowSizes;
            insertQrow += gappedExtension_rowSizes;
        }

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
            previousInsertQ = *insertQrow;

			subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;

            // Check for scoring row wrap-around
            if (matchRow < gappedExtension_matchRow)
            {
                matchRow += gappedExtension_rowSizes;
                insertSrow += gappedExtension_rowSizes;
                insertQrow += gappedExtension_rowSizes;
            }
		}

        // -----CELLS LEFT OF ROW DROPOFF -----
		if (!(dropoffThreshold > previousMatch &&
			dropoffThreshold > previousInsertS &&
			dropoffThreshold > previousInsertQ))
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

                // Check for scoring row wrap-around
                if (matchRow < gappedExtension_matchRow)
                {
                    matchRow += gappedExtension_rowSizes;
                    insertSrow += gappedExtension_rowSizes;
                    insertQrow += gappedExtension_rowSizes;
                }

				subjectDistance++;
			}
		}

		// Record dropoff position
		rowDropoff = subjectPosition + 1;

//		print(gappedExtension_matchRow, subject, rowDropoff, columnDropoff);

		queryDistance++;
	}

	dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = gappedExtension_traceback;

	return dpResults;
}

// Perform dynamic programming to explore possible END points and alignments that start at
// the given seed
struct dpResults gappedExtension_dpAfterSeed(struct PSSMatrix PSSMatrix, int4 dropoff,
       struct unpackRegion* unpackRegion, int4 subjectLength, int4 seedSubjectOffset)
{
	int2 **queryPosition, **bestQueryPosition, **queryEnd;
	int2 *matrixColumn;
	unsigned char *subject, *rowDropoff, *columnDropoff, *maxRowDropoff, *newSubject;
	unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
	unsigned char **tracebackRow, *tracebackColumn;
	int4 bestScore = 0, dropoffThreshold;
	int4 *matchRow, *insertQrow, *insertSrow, rowOffset, *endMatchRow;
	int4 queryDistance, subjectDistance;
	int4 oldMatch, match, previousOldMatch, previousOldInsertS, previousOldInsertQ;
	int4 previousMatch, previousInsertS, previousInsertQ = 0;
	struct dpResults dpResults;
	unsigned char leftOfDropoff;
	int4 queryLength;
    unsigned char *tracebackDataPtr, *newTracebackData;

    subject = unpackRegion->unpackedSubject + seedSubjectOffset;

    queryLength = PSSMatrix.length;
	subjectEnd = subject + subjectLength;
	queryEnd = PSSMatrix.matrix + queryLength;
    endMatchRow = gappedExtension_matchRow + gappedExtension_rowSizes;

	// Determine lowest score before dropoff
	dropoffThreshold = -dropoff;

    // If more rows are required
    if (queryLength > gappedExtension_numRows)
    {
        // Increase number of row pointers
    	gappedExtension_numRows = queryLength;
        gappedExtension_traceback = (unsigned char**)global_realloc(gappedExtension_traceback,
                                    sizeof(unsigned char*) * (gappedExtension_numRows));
    }

    tracebackDataPtr = gappedExtension_tracebackData;

    bestSubjectPosition = subjectPosition = subject + 1;
	bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

	// Initialize rows
	matchRow = gappedExtension_matchRow + 1;
	insertQrow = gappedExtension_insertQrow + 1;
	insertSrow = gappedExtension_insertSrow + 1;

	// Set initial row dropoff and column dropoff
	rowDropoff = subject + subjectLength - 1;
	columnDropoff = subject;

    // Calculate maximum row dropoff after this row is processed
    maxRowDropoff = columnDropoff + ((dropoff / parameters_extendGap) + 2);
    if (maxRowDropoff > subjectEnd)
        maxRowDropoff = subjectEnd;

    // Initialize traceback pointer for this row
	tracebackRow = gappedExtension_traceback + (queryPosition - PSSMatrix.matrix);
	*tracebackRow = tracebackDataPtr - (subjectPosition - subject);
    tracebackColumn = tracebackDataPtr;
    tracebackDataPtr += (maxRowDropoff - subjectPosition + 1);

    // Unpack more of the subject if required
    unpack_extendRegionEnd(maxRowDropoff - subject + seedSubjectOffset, unpackRegion);

    // If the unpacked subject has moved in memory
    if (subject != unpackRegion->unpackedSubject + seedSubjectOffset)
    {
        // Update all pointers to point at the new subject
        newSubject = unpackRegion->unpackedSubject + seedSubjectOffset;
        subjectPosition = (subjectPosition - subject) + newSubject;
        rowDropoff = (rowDropoff - subject) + newSubject;
        columnDropoff = (columnDropoff - subject) + newSubject;
        bestSubjectPosition = (bestSubjectPosition - subject) + newSubject;
        subject = newSubject;
        subjectEnd = subject + subjectLength;
    }

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
//    print2(gappedExtension_matchRow, subject, rowDropoff, columnDropoff);

	queryDistance = 0;
	queryPosition++; tracebackRow++;

	// -----REMAINING ROWS-----
	while (queryPosition < queryEnd && rowDropoff > columnDropoff)
	{
		subjectPosition = columnDropoff + 1;

        // Calculate maximum row dropoff after this row is processed
        maxRowDropoff = rowDropoff + ((PSSMatrix.highestValue / parameters_extendGap) + 2);
        if (maxRowDropoff > subjectEnd)
        	maxRowDropoff = subjectEnd;

        // If not enough space in traceback data for this row, realloc
        if (maxRowDropoff - subjectPosition >= gappedExtension_tracebackAlloc -
            (tracebackDataPtr - gappedExtension_tracebackData))
        {
			newTracebackData = (unsigned char*)global_malloc(sizeof(char) * gappedExtension_tracebackAlloc * 2);

            // Move existing data to end of new block
            memcpy(newTracebackData, gappedExtension_tracebackData, gappedExtension_tracebackAlloc);

            // Update traceback data pointer
			tracebackDataPtr = newTracebackData + (tracebackDataPtr - gappedExtension_tracebackData);

            // Update existing rows to point to new data block
			while (tracebackRow > gappedExtension_traceback)
            {
				tracebackRow--;
            	*tracebackRow = newTracebackData + (*tracebackRow - gappedExtension_tracebackData);
            }

            // Data block is now double the size
            gappedExtension_tracebackAlloc *= 2;

            free(gappedExtension_tracebackData);
            gappedExtension_tracebackData = newTracebackData;
        }

        // Initialize traceback pointer for this row
        tracebackRow = gappedExtension_traceback + (queryPosition - PSSMatrix.matrix);
        *tracebackRow = tracebackDataPtr - (subjectPosition - subject);
        tracebackColumn = tracebackDataPtr;
        tracebackDataPtr += (maxRowDropoff - subjectPosition + 1);

        // Unpack more of the subject if required
        unpack_extendRegionEnd(maxRowDropoff - subject + seedSubjectOffset, unpackRegion);

        // If the unpacked subject has moved in memory
        if (subject != unpackRegion->unpackedSubject + seedSubjectOffset)
        {
        	// Update all pointers to point at the new subject
            newSubject = unpackRegion->unpackedSubject + seedSubjectOffset;
        	subjectPosition = (subjectPosition - subject) + newSubject;
			rowDropoff = (rowDropoff - subject) + newSubject;
			columnDropoff = (columnDropoff - subject) + newSubject;
            bestSubjectPosition = (bestSubjectPosition - subject) + newSubject;
            subject = newSubject;
        	subjectEnd = subject + subjectLength;
        }

        // Reset rows
		rowOffset = (subjectPosition - subject);
        matchRow = gappedExtension_matchRow + (rowOffset % gappedExtension_rowSizes);
        insertQrow = gappedExtension_insertQrow + (rowOffset % gappedExtension_rowSizes);
        insertSrow = gappedExtension_insertSrow + (rowOffset % gappedExtension_rowSizes);

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

        // Check for scoring rows wrap-around
        if (matchRow >= endMatchRow)
        {
            matchRow -= gappedExtension_rowSizes;
            insertQrow -= gappedExtension_rowSizes;
            insertSrow -= gappedExtension_rowSizes;
        }

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
            previousInsertQ = *insertQrow;

			subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;

            // Check for scoring rows wrap-around
            if (matchRow >= endMatchRow)
            {
            	matchRow -= gappedExtension_rowSizes;
            	insertQrow -= gappedExtension_rowSizes;
            	insertSrow -= gappedExtension_rowSizes;
            }
		}

		// -----CELLS RIGHT OF ROW DROPOFF -----
		if (!(dropoffThreshold > previousMatch &&
			dropoffThreshold > previousInsertS &&
			dropoffThreshold > previousInsertQ))
		{
			while (subjectPosition < subjectEnd)
			{
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

                // Check for scoring rows wrap-around
                if (matchRow >= endMatchRow)
                {
                    matchRow -= gappedExtension_rowSizes;
                    insertQrow -= gappedExtension_rowSizes;
                    insertSrow -= gappedExtension_rowSizes;
                }

				subjectDistance++;
			}
		}

		// Record dropoff position
		rowDropoff = subjectPosition - 1;

		queryDistance++;
		queryPosition++; tracebackRow++;

//        if (dloc==88197331)
//			gappedExtension_printAfterRow(gappedExtension_matchRow, subject, rowDropoff,
//                                          columnDropoff);
	}

	dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = gappedExtension_traceback;

	return dpResults;
}

void gappedExtension_free()
{
	// Free memory used by traceback array
    free(gappedExtension_traceback);
    free(gappedExtension_tracebackData);

//    printf("gappedExtension_tracebackAlloc=%d\n", gappedExtension_tracebackAlloc);

    // Free memory used to store row scores
    free(gappedExtension_matchRow);
    free(gappedExtension_insertQrow);
    free(gappedExtension_insertSrow);
}

// Debugging routine
void gappedExtension_printBeforeRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                    unsigned char* columnDropoff)
{
	// TODO: Adjust to work with wrap-around
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
void gappedExtension_printAfterRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                   unsigned char* columnDropoff)
{
	// TODO: Adjust to work with wrap-around
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
