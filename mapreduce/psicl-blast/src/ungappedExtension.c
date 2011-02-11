// ungappedExtension.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code to perform an ungapped extension

#include "blast.h"

unsigned char* ungappedExtension_subjectEndReached;
uint4 ungappedExtension_subjectEndReachedFP;
int4 ungappedExtension_bestScore;
struct memBlocks* ungappedExtension_extensions;

int4 ungappedExtension_minus3reward;
int4 ungappedExtension_tableMatchesReward;

// Prototypes
struct coordinate ungappedExtension_findNucleotideSeed(struct ungappedExtension* ungappedExtension,
                                                       struct PSSMatrix PSSMatrix, unsigned char* subject);
// Shucai
struct coordinate ungappedExtension_findProteinSeed(struct ungappedExtension* ungappedExtension,
       struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, unsigned char* subject);

// Initialize the creation of ungapped extensions
void ungappedExtension_initialize()
{
	ungappedExtension_extensions = memBlocks_initialize(sizeof(struct ungappedExtension),
                                   constants_initialAllocUngappedExtensions);
	ungappedExtension_minus3reward = parameters_matchScore * -3;
    ungappedExtension_tableMatchesReward = parameters_matchScore * parameters_wordTableLetters;
}

// Perform an ungapped extension between points queryStart,subjectStart and queryEnd,subjectEnd
// and extend in each direction until score drops below best score yet minus a dropoff parameter
// Shucai
struct ungappedExtension* ungappedExtension_extend(int4 queryoffset, unsigned char* subjectHit,
	uint4 lastHitFP, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, 
	unsigned char* subject, unsigned char *startAddressFP)
{
	//Shucai
	int2 *queryPosition;
	unsigned char *subjectPosition, *subjectStart, *subjectEnd;
	int4 changeSinceBest = 0;
	int4 dropoff, originalDropoff;

    originalDropoff = dropoff = -statistics_ungappedNominalDropoff;
	ungappedExtension_bestScore = 0;

	// Start at queryEnd,subjectEnd (right/last hit position)
	queryPosition = PSSMatrixFP.matrix + queryoffset * encoding_numCodes;
	subjectPosition = subjectStart = subjectHit;

	// Extend the start of the hit backwards until dropoff
	while (changeSinceBest > dropoff)
	{
		//changeSinceBest += (*queryPosition)[*subjectPosition];
		changeSinceBest += queryPosition[*subjectPosition];

        // If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				//Shucai
				queryPosition -= encoding_numCodes; 
				subjectPosition--;
				//Shucai
				changeSinceBest = queryPosition[*subjectPosition];
			}
			while (changeSinceBest > 0);

			subjectStart = subjectPosition;
		}
		//Shucai
		queryPosition -= encoding_numCodes; 
		subjectPosition--;
	}

	// Correct for extra decrement
	subjectStart++;

	// If best start point is right of previous hit which helped trigger this extension
	// then stop now
	// Shucai
	//if (subjectStart - startAddressFP > lastHitFP)
	if (subjectStart -  subject > lastHitFP)
	{
		//Shucai
		//ungappedExtension_subjectEndReachedFP = subjectHit - startAddressFP;
		ungappedExtension_subjectEndReachedFP = subjectHit - subject;
		return NULL;
	}

	// Starting at right/last hit position again
	//Shucai
	queryPosition = PSSMatrixFP.matrix + (queryoffset + 1) * encoding_numCodes;
    subjectEnd = subjectHit;

	subjectPosition = subjectHit + 1;
    changeSinceBest = 0;

    // May need to alter dropoff so we also dropoff if below zero
    if (-ungappedExtension_bestScore > originalDropoff)
    {
    	dropoff = -ungappedExtension_bestScore;
    }

	// Extend end of alignment until dropoff
	while (changeSinceBest > dropoff)
	{
		//Shucai
		changeSinceBest += queryPosition[*subjectPosition];

        // If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				//Shucai
				queryPosition += encoding_numCodes; 
				subjectPosition++;
				//Shucai
				changeSinceBest = queryPosition[*subjectPosition];
			}
			while (changeSinceBest > 0);

			subjectEnd = subjectPosition;

			// Check need for change in dropoff
            if ((dropoff = -ungappedExtension_bestScore) < originalDropoff)
            {
            	dropoff = originalDropoff;
            }
        }
		//Shucai
		queryPosition += encoding_numCodes; 
		subjectPosition++;
	}

	// Correct for extra increment
	subjectEnd--;
	//Shucai
	//ungappedExtension_subjectEndReachedFP = subjectEnd - startAddressFP;
	ungappedExtension_subjectEndReachedFP = subjectEnd - subject;

    // If extension scored above trigger for gapping, create object and return it
    if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger)
    {
    	int4 diagonal;
        struct ungappedExtension* newUngappedExtension;
        newUngappedExtension = memBlocks_newEntry(ungappedExtension_extensions);

        // Calculate diagonal
        // Shucai
		diagonal = (subjectHit - subject) - queryoffset;

        // Determine offsets from pointers
        newUngappedExtension->start.subjectOffset = subjectStart - subject;
        newUngappedExtension->end.subjectOffset = subjectEnd - subject;
        newUngappedExtension->start.queryOffset = newUngappedExtension->start.subjectOffset - diagonal;
        newUngappedExtension->end.queryOffset = newUngappedExtension->end.subjectOffset - diagonal;

        // Find the seed point
        newUngappedExtension->seed = ungappedExtension_findProteinSeed(newUngappedExtension, PSSMatrix, PSSMatrixFP, subject);
        // Initialize next to null
        newUngappedExtension->next = NULL;
        newUngappedExtension->nominalScore = ungappedExtension_bestScore;
        newUngappedExtension->status = ungappedExtension_UNGAPPED;

        return newUngappedExtension;
    }
    else
    {
    	return NULL;
    }
}

// Perform an ungapped extension when the seed is only a single hit on the diagonal, rather
// than a pair of hits.
// Shucai
struct ungappedExtension* ungappedExtension_oneHitExtend(int4 queryoffset,
	unsigned char* subjectHit, struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, 
	unsigned char* subject, unsigned char *startAddressFP)
{
	//Shucai
	int2* queryPosition;
	unsigned char* subjectPosition, *subjectStart, *subjectEnd;
	int4 changeSinceBest = 0;
	int4 dropoff, originalDropoff;

	originalDropoff = dropoff = -statistics_ungappedNominalDropoff;
    ungappedExtension_bestScore = 0;

	// Start at queryEnd,subjectEnd (right/last hit position)
	// Shucai
	queryPosition = PSSMatrixFP.matrix + queryoffset * encoding_numCodes;
	
	subjectPosition = subjectStart = subjectHit;

	// Extend the start of the hit forwards until dropoff
	while (changeSinceBest > dropoff)
	{
		//Shucai
		changeSinceBest += queryPosition[*subjectPosition];
		
		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				//Shucai
				queryPosition = queryPosition - encoding_numCodes; 
				subjectPosition--;
				//Shucai
				changeSinceBest = queryPosition[*subjectPosition];
			}
			while (changeSinceBest > 0);

			subjectStart = subjectPosition;
		}
		//Shucai
		queryPosition = queryPosition - encoding_numCodes; 
		subjectPosition--;
	}

	// Correct for extra decrement
	subjectStart++;

	// Starting at right/last hit position again
	//Shucai
	queryPosition = PSSMatrixFP.matrix + (queryoffset + 1) * encoding_numCodes;
	subjectPosition = subjectEnd = subjectHit + 1;
	changeSinceBest = 0;

    // May need to alter dropoff so we also dropoff if below zero
    if (-ungappedExtension_bestScore > originalDropoff)
    {
    	dropoff = -ungappedExtension_bestScore;
    }

	// Extend end of alignment until dropoff
	while (changeSinceBest > dropoff)
	{
		//Shucai
		changeSinceBest += queryPosition[*subjectPosition];
		
		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				//Shucai
				queryPosition = queryPosition + encoding_numCodes; 
				subjectPosition++;
				//Shucai
				changeSinceBest = queryPosition[*subjectPosition];
			}
			while (changeSinceBest > 0);

			subjectEnd = subjectPosition;

			// Check need for change in dropoff
            if ((dropoff = -ungappedExtension_bestScore) < originalDropoff)
            {
            	dropoff = originalDropoff;
            }
        }
		//Shucai
		queryPosition = queryPosition + encoding_numCodes; 
		subjectPosition++;
	}

	// Correct for extra increment
	subjectEnd--;

    // Record the point we got to extending forwards
	// Shucai
	//ungappedExtension_subjectEndReachedFP = subjectPosition - startAddressFP;
	ungappedExtension_subjectEndReachedFP = subjectPosition - subject;

    // If extension scored above trigger for gapping, create object and return it
    if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger)
    {
    	int4 diagonal;
        struct ungappedExtension* newUngappedExtension;
        newUngappedExtension = memBlocks_newEntry(ungappedExtension_extensions);

        // Calculate diagonal
		//Shucai
		diagonal = (subjectHit - subject) - queryoffset;

        // Determine offsets from pointers
        newUngappedExtension->start.subjectOffset = subjectStart - subject;
        newUngappedExtension->end.subjectOffset = subjectEnd - subject;
        newUngappedExtension->start.queryOffset = newUngappedExtension->start.subjectOffset - diagonal;
        newUngappedExtension->end.queryOffset = newUngappedExtension->end.subjectOffset - diagonal;

        // Find the seed point
        newUngappedExtension->seed = ungappedExtension_findProteinSeed(newUngappedExtension, PSSMatrix, PSSMatrixFP, subject);
        // Initialize next to null
        newUngappedExtension->next = NULL;
        newUngappedExtension->nominalScore = ungappedExtension_bestScore;
        newUngappedExtension->status = ungappedExtension_UNGAPPED;

        return newUngappedExtension;
    }
    else
    {
    	return NULL;
    }
}

// Perform one-hit seeded ungapped extension for nucleotide, 1 packed-byte at a time
struct ungappedExtension* ungappedExtension_nucleotideExtend(int4 queryHitOffset,
	int4 subjectHitOffset, struct PSSMatrix PSSMatrix, unsigned char* subject,
    uint4 subjectLength)
{
	unsigned char* queryPosition, *minQueryPosition, *maxQueryPosition;
	unsigned char* subjectPosition, *subjectStart, *subjectEnd;
	int4 dropoff, originalDropoff;
    int4 changeSinceBest = 0;
    int4 matchLettersScore;

	originalDropoff = dropoff = -statistics_ungappedNominalDropoff;

    // Start with score for lookup-table nucleotide match that is not aligned
    ungappedExtension_bestScore = ungappedExtension_tableMatchesReward;

    // Determine minimum query position; either start of the query or start of the second strand
    if (queryHitOffset <= PSSMatrix.strandLength)
    {
        if (queryHitOffset < subjectHitOffset * 4)
            minQueryPosition = PSSMatrix.bytePackedCodes;
        else
            minQueryPosition = PSSMatrix.bytePackedCodes + queryHitOffset - subjectHitOffset * 4;
	}
    else
    {
        if (queryHitOffset - PSSMatrix.strandLength < subjectHitOffset * 4)
            minQueryPosition = PSSMatrix.bytePackedCodes + PSSMatrix.strandLength;
        else
            minQueryPosition = PSSMatrix.bytePackedCodes + queryHitOffset - subjectHitOffset * 4;
    }

	// Start left of hit location
	queryPosition = PSSMatrix.bytePackedCodes + queryHitOffset - parameters_wordTableLetters - 4;
	subjectPosition = subjectStart = subject + subjectHitOffset - parameters_wordTableBytes - 1;

    // Consider partial match of first byte before hit
	matchLettersScore = PSSMatrix_packedLeftMatchScores[*queryPosition ^ *subjectPosition];
    ungappedExtension_bestScore += matchLettersScore;
	changeSinceBest = -matchLettersScore;

    // Move back through alignment until start of query or subject, or until dropoff
    while (queryPosition > minQueryPosition)
    {
    	// Add score of matching entire bytes
		changeSinceBest += PSSMatrix_packedScore[*queryPosition ^ *subjectPosition];

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
        {
        	printf("<%d< ", PSSMatrix_packedScore[*queryPosition ^ *subjectPosition]);
            printf("["); encoding_printLetters(*queryPosition, 4);
            printf(","); encoding_printLetters(*subjectPosition, 4); printf("]\n");
		}
        #endif

        // If we possibly have a new best score
        if (changeSinceBest > ungappedExtension_minus3reward)
        {
            // Get score for matching individual letters in next byte
        	queryPosition-=4; subjectPosition--;
	        matchLettersScore = PSSMatrix_packedLeftMatchScores[*queryPosition ^ *subjectPosition];

            // If best score
            if (changeSinceBest + matchLettersScore > 0)
            {
                // Mark new best position
                subjectStart = subjectPosition;

                // Update best score and change since best
                ungappedExtension_bestScore += changeSinceBest + matchLettersScore;
                changeSinceBest = -matchLettersScore;

                #ifdef VERBOSE
                if (parameters_verboseDloc == blast_dloc)
                    printf("(Best=%d)\n", ungappedExtension_bestScore);
                #endif
            }
        }
        else
        {
        	// Decrease in score, check dropoff
			if (changeSinceBest < dropoff)
            	break;

            queryPosition-=4; subjectPosition--;
        }
    }

    // Determine maximum query position; either end of the query or end of the first strand
    if (queryHitOffset <= PSSMatrix.strandLength)
    {
        if (PSSMatrix.strandLength - queryHitOffset < subjectLength - subjectHitOffset * 4)
            maxQueryPosition = PSSMatrix.bytePackedCodes + PSSMatrix.strandLength - 4;
        else
            maxQueryPosition = PSSMatrix.bytePackedCodes + (subjectLength - subjectHitOffset * 4)
                             + queryHitOffset - 4;
	}
    else
    {
        if (PSSMatrix.length - queryHitOffset < subjectLength - subjectHitOffset * 4)
            maxQueryPosition = PSSMatrix.bytePackedCodes + PSSMatrix.length - 4;
        else
            maxQueryPosition = PSSMatrix.bytePackedCodes + (subjectLength - subjectHitOffset * 4)
                             + queryHitOffset - 4;
    }

    // Starting right of hit position
	queryPosition = PSSMatrix.bytePackedCodes + queryHitOffset;
	subjectPosition = subjectEnd = subject + subjectHitOffset;
	changeSinceBest = 0;

    // May need to alter dropoff so we also dropoff if below zero
    if (-ungappedExtension_bestScore > originalDropoff)
    {
    	dropoff = -ungappedExtension_bestScore;
    }

    // Consider partial match of first byte after hit
	matchLettersScore = PSSMatrix_packedRightMatchScores[*queryPosition ^ *subjectPosition];
    ungappedExtension_bestScore += matchLettersScore;
	changeSinceBest = -matchLettersScore;

    // Move forward through alignment until end of query or subject, or until dropoff
    while (queryPosition < maxQueryPosition)
    {
		// Score of matching entire bytes
		changeSinceBest += PSSMatrix_packedScore[*queryPosition ^ *subjectPosition];

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
        {
        	printf(">%d> ", PSSMatrix_packedScore[*queryPosition ^ *subjectPosition]);
            printf("["); encoding_printLetters(*queryPosition, 4);
            printf(","); encoding_printLetters(*subjectPosition, 4); printf("]\n");
            printf("changeSinceBest=%d\n", changeSinceBest);
		}
        #endif

        // If we possibly have a new best score
        if (changeSinceBest > ungappedExtension_minus3reward)
        {
            // Get score for matching individual letters in next byte
        	queryPosition+=4; subjectPosition++;
	        matchLettersScore = PSSMatrix_packedRightMatchScores[*queryPosition ^ *subjectPosition];

            // If best score
            if (changeSinceBest + matchLettersScore > 0)
            {
                // Mark new best position
                subjectEnd = subjectPosition;

                // Update best score and change since best
                ungappedExtension_bestScore += changeSinceBest + matchLettersScore;
                changeSinceBest = -matchLettersScore;

                #ifdef VERBOSE
                if (parameters_verboseDloc == blast_dloc)
                    printf("(Best=%d)\n", ungappedExtension_bestScore);
                #endif
            }
        }
        else
        {
        	// Decrease in score, check dropoff
			if (changeSinceBest < dropoff)
            	break;

            queryPosition+=4; subjectPosition++;
        }
    }

    // Record the point we got to extending forwards
    ungappedExtension_subjectEndReached = subjectPosition;

    // If extension scored above trigger for gapping, create object and return it
    if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger)
    {
    	int4 diagonal;
        struct ungappedExtension* newUngappedExtension;
        newUngappedExtension = memBlocks_newEntry(ungappedExtension_extensions);

        // Correct for extra decrement
        subjectStart++;
        // Correct for extra increment
        subjectEnd--;

        // Calculate diagonal
        diagonal = subjectHitOffset * 4 - queryHitOffset;

        // Determine offsets from pointers
        newUngappedExtension->start.subjectOffset = (subjectStart - subject) * 4;
        newUngappedExtension->end.subjectOffset = (subjectEnd - subject) * 4;
        newUngappedExtension->start.queryOffset = newUngappedExtension->start.subjectOffset - diagonal;
        newUngappedExtension->end.queryOffset = newUngappedExtension->end.subjectOffset - diagonal;

		newUngappedExtension->seed.queryOffset = -1;
		newUngappedExtension->seed.subjectOffset = -1;

        // Initialize next to null
        newUngappedExtension->next = NULL;
        newUngappedExtension->nominalScore = ungappedExtension_bestScore;
        newUngappedExtension->status = ungappedExtension_UNGAPPED;

        #ifdef VERBOSE
        if (parameters_verboseDloc == blast_dloc)
        {
            printf("Hit=%d,%d\n", queryHitOffset, subjectHitOffset);
            printf("%d,%d - %d,%d\n", newUngappedExtension->start.queryOffset, newUngappedExtension->start.subjectOffset,
                                      newUngappedExtension->end.queryOffset, newUngappedExtension->end.subjectOffset);
                                      fflush(stdout);
            printf("seed=%d,%d\n", newUngappedExtension->seed.queryOffset, newUngappedExtension->seed.subjectOffset);
		}
		#endif

        return newUngappedExtension;
    }
    else
    {
    	return NULL;
    }
}

// Check that the given coordinate pair is part of a hit of length at least hitLength
int ungappedExtension_checkHit(uint4 queryHitOffset, uint4 subjectHitOffset,
                               struct PSSMatrix PSSMatrix, unsigned char* subject,
                               uint4 subjectLength, uint4 hitLength)
{
	uint4 match, match4, targetScore, stopScore;
	unsigned char* queryPosition, *minQueryPosition, *maxQueryPosition;
	unsigned char* subjectPosition;

    ungappedExtension_bestScore = 0;
    match4 = parameters_matchScore * 4;
    stopScore = (hitLength - 4) * parameters_matchScore;
    targetScore = hitLength * parameters_matchScore;

    // Determine minimum query position; either start of the query or start of the second strand
    if (queryHitOffset < PSSMatrix.strandLength)
    {
        if (queryHitOffset < subjectHitOffset * 4)
            minQueryPosition = PSSMatrix.bytePackedCodes;
        else
            minQueryPosition = PSSMatrix.bytePackedCodes + queryHitOffset - subjectHitOffset * 4;
	}
    else
    {
        if (queryHitOffset - PSSMatrix.strandLength < subjectHitOffset * 4)
            minQueryPosition = PSSMatrix.bytePackedCodes + PSSMatrix.strandLength;
        else
            minQueryPosition = PSSMatrix.bytePackedCodes + queryHitOffset - subjectHitOffset * 4;
    }

	// Start left of hit location
	queryPosition = PSSMatrix.bytePackedCodes + queryHitOffset;
	subjectPosition = subject + subjectHitOffset;

    while (queryPosition >= minQueryPosition)
    {
    	// Add score of matching next byte
		match = PSSMatrix_packedLeftMatchScores[*queryPosition ^ *subjectPosition];

//        printf("<%d< (%d)", match, queryPosition - PSSMatrix.bytePackedCodes);

        if (match == match4 && ungappedExtension_bestScore < stopScore)
        {
			ungappedExtension_bestScore += match;
		}
        else
        {
        	if (match > 0)
            	ungappedExtension_bestScore += match;
        	break;
        }

        queryPosition-=4; subjectPosition--;
    }

    if (ungappedExtension_bestScore >= targetScore)
    	return ungappedExtension_bestScore;

    // Determine maximum query position; either end of the query or end of the first strand
    if (queryHitOffset < PSSMatrix.strandLength)
    {
        if (PSSMatrix.strandLength - queryHitOffset < subjectLength - subjectHitOffset * 4)
            maxQueryPosition = PSSMatrix.bytePackedCodes + PSSMatrix.strandLength - 4;
        else
            maxQueryPosition = PSSMatrix.bytePackedCodes + (subjectLength - subjectHitOffset * 4)
                             + queryHitOffset - 4;
	}
    else
    {
        if (PSSMatrix.length - queryHitOffset < subjectLength - subjectHitOffset * 4)
            maxQueryPosition = PSSMatrix.bytePackedCodes + PSSMatrix.length - 4;
        else
            maxQueryPosition = PSSMatrix.bytePackedCodes + (subjectLength - subjectHitOffset * 4)
                             + queryHitOffset - 4;
    }

    // Starting right of hit position
	queryPosition = PSSMatrix.bytePackedCodes + queryHitOffset + 4;
	subjectPosition = subject + subjectHitOffset + 1;

    while (queryPosition <= maxQueryPosition)
    {
    	// Add score of matching next byte
		match = PSSMatrix_packedRightMatchScores[*queryPosition ^ *subjectPosition];

//        printf(">%d> (%d)", match, queryPosition - PSSMatrix.bytePackedCodes);

        if (match == match4 && ungappedExtension_bestScore < stopScore)
        {
			ungappedExtension_bestScore += match;
		}
        else
        {
        	if (match > 0)
            	ungappedExtension_bestScore += match;
        	break;
        }

        queryPosition+=4; subjectPosition++;
    }

    // Record the point we got to extending forwards
    ungappedExtension_subjectEndReached = subjectPosition;

    if (ungappedExtension_bestScore >= targetScore)
    	return ungappedExtension_bestScore;
	else
    	return 0;
}

// Find seed point for an ungapped extension
void ungappedExtension_findSeed(struct ungappedExtension* ungappedExtension,
	struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, unsigned char* subject)
{
	// Find seed if not already recorded
	if (ungappedExtension->seed.queryOffset == -1 && ungappedExtension->seed.subjectOffset == -1)
    {
        if (encoding_alphabetType == encoding_protein)
        {
            ungappedExtension->seed
            	= ungappedExtension_findProteinSeed(ungappedExtension, PSSMatrix, PSSMatrixFP, subject);
        }
        else
        {
            ungappedExtension->seed
            	= ungappedExtension_findNucleotideSeed(ungappedExtension, PSSMatrix, subject);
        }
	}
}

// Find the seed point (middle of the longest run of matching bytes) of an ungapped extension
// for a bytepacked nucleotide sequence
struct coordinate ungappedExtension_findNucleotideSeed(struct ungappedExtension* ungappedExtension,
                                                       struct PSSMatrix PSSMatrix, unsigned char* subject)
{
	unsigned char *queryPosition, *subjectPosition, *subjectEnd;
    unsigned char *queryStartRun, *subjectStartRun;
    int4 queryMiddle, subjectMiddle;
    int4 run = 0, longestRun = 0, diagonal;
	struct coordinate seed;

    // Start at beginning of ungapped alignment
    queryPosition = PSSMatrix.bytePackedCodes + ungappedExtension->start.queryOffset;
    subjectPosition = subject + ungappedExtension->start.subjectOffset / 4;
	queryStartRun = queryPosition;
	subjectStartRun = subjectPosition;
    subjectMiddle = (ungappedExtension->start.subjectOffset + ungappedExtension->end.subjectOffset) / 8;

    subjectEnd = subject + ungappedExtension->end.subjectOffset / 4;
    // Slide across ungapped alignment until end
    while (subjectPosition <= subjectEnd)
    {
		if (*subjectPosition == *queryPosition)
        {
	    	// Bytes match, increase length of run
        	run++;
        }
        else
        {
        	// Bytes mismatch, end of run, check if the longest run
        	if (run > longestRun)
            {
            	// If longest run find its middle
				longestRun = run;
				subjectMiddle = ((subjectPosition - 1 - subject) + (subjectStartRun - subject)) / 2;
            }
            // Reset length and start of run
        	run = 0;
            queryStartRun = queryPosition;
            subjectStartRun = subjectPosition;
        }

        subjectPosition++;
        queryPosition+=4;
    }

    // Calculate the query middle based on the subject middle
    diagonal = ungappedExtension->start.subjectOffset - ungappedExtension->start.queryOffset;
    subjectMiddle *= 4;
    queryMiddle = subjectMiddle - diagonal;

    // Use the middle of the seed byte
    seed.queryOffset = queryMiddle + 2;
    seed.subjectOffset = subjectMiddle + 2;

    return seed;
}

// Find the seed point (the middle of the 11-aa long highest scoring region) of the
// given ungapped extension
struct coordinate ungappedExtension_findProteinSeed(struct ungappedExtension* ungappedExtension,
       struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, unsigned char* subject)
{
	//Shucai
	//int2 **queryWindowStart, **queryWindowEnd;
	int2 *queryWindowStart, *queryWindowEnd;
	unsigned char *subjectWindowStart, *subjectWindowEnd;
	//Shucai
	//int2** bestQueryPosition;
	int2* bestQueryPosition;
	unsigned char* bestSubjectPosition;
	int4 bestSegmentScore;
	int4 nominalScore, count;
	struct coordinate seed;

    // If the length of the ungapped extension is 11 or less
    if (ungappedExtension->end.queryOffset - ungappedExtension->start.queryOffset < 11)
    {
        // The seed point is the middle of the extension
        seed.queryOffset = (ungappedExtension->end.queryOffset +
                            ungappedExtension->start.queryOffset) / 2;
        seed.subjectOffset = (ungappedExtension->end.subjectOffset +
                              ungappedExtension->start.subjectOffset) / 2;
    }
    else
    {
        // Else find the highest scoring length-11 segment of the ungapped extension
        //Shucai
		//queryWindowStart = queryWindowEnd = PSSMatrix.matrix + ungappedExtension->start.queryOffset;
		queryWindowStart = queryWindowEnd = PSSMatrixFP.matrix + ungappedExtension->start.queryOffset * encoding_numCodes;
        subjectWindowStart = subjectWindowEnd = subject + ungappedExtension->start.subjectOffset;

        // Find initial score for first 11 positions
        nominalScore = 0;
        count = 0;
        while (count < 11)
        {
            // Add to tally score for query/subject at this position
			// Shucai
			//nominalScore += (*queryWindowEnd)[*subjectWindowEnd];
            nominalScore += queryWindowEnd[*subjectWindowEnd];
			
			//Shucai
            //queryWindowEnd++;
			queryWindowEnd += encoding_numCodes;
            subjectWindowEnd++;
            count++;
        }
		
		//Shucai
        //queryWindowEnd--;
		queryWindowEnd -= encoding_numCodes;
        subjectWindowEnd--;

        // By default first-11 positions gives best position and score
        bestQueryPosition = queryWindowStart;
        bestSubjectPosition = subjectWindowStart;
        bestSegmentScore = nominalScore;

        // Now slide the window across and record the better scores/positions
        // Shucai
		//while (queryWindowEnd < PSSMatrix.matrix + ungappedExtension->end.queryOffset)
		while (queryWindowEnd < PSSMatrixFP.matrix + ungappedExtension->end.queryOffset * encoding_numCodes)
        {
            // Advance window end, add new position value
			// Shucai
			//queryWindowEnd++;
            queryWindowEnd += encoding_numCodes;
            subjectWindowEnd++;
			//Shucai
			//nominalScore += (*queryWindowEnd)[*subjectWindowEnd];
            nominalScore += queryWindowEnd[*subjectWindowEnd];

            // Remove position that we will leave behind
            // Shucai
			//nominalScore -= (*queryWindowStart)[*subjectWindowStart];
			nominalScore -= queryWindowStart[*subjectWindowStart];

            // Advance window start
            // Shucai
			//queryWindowStart++;
			queryWindowStart += encoding_numCodes;
            subjectWindowStart++;

            // Check if best window position yet
            if (nominalScore > bestSegmentScore)
            {
                bestSegmentScore = nominalScore;
                bestQueryPosition = queryWindowStart;
                bestSubjectPosition = subjectWindowStart;
            }
        }

        // Middle of the best window is the seed position
        //Shucai
		//seed.queryOffset = (bestQueryPosition - PSSMatrix.matrix) + 5;
		seed.queryOffset = (bestQueryPosition - PSSMatrixFP.matrix) / encoding_numCodes + 5;
        seed.subjectOffset = bestSubjectPosition + 5 - subject;
    }

	return seed;
}

// Print details of an ungapped extension (for debugging purposes)
void ungappedExtension_print(struct ungappedExtension* extension)
{
	printf("Ungapped extension %4d,%4d to %4d,%4d score=%4d status=%d seed=%4d,%4d\n",
	extension->start.queryOffset, extension->start.subjectOffset,
	extension->end.queryOffset, extension->end.subjectOffset,
	extension->nominalScore, extension->status,
    extension->seed.queryOffset, extension->seed.subjectOffset);
}


struct ungappedExtension* ungappedExtension_simpleScoring(
       unsigned char* sequence1, uint4 offset1, unsigned char* sequence2, uint4 offset2)
{
	int4 position1, position2, start;
	int4 dropoff, score;

	dropoff = 5;
    ungappedExtension_bestScore = 0;
    score = 0;

	position1 = offset1;
	position2 = start = offset2;

    while (position2 > 0 && position1 > 0 && score + dropoff >= ungappedExtension_bestScore)
    {
    	if (sequence1[position1] == sequence2[position2])
        {
            score++; printf("1");
		}
        else
		{
            score--; printf("0");
		}

        if (score > ungappedExtension_bestScore)
        {
			ungappedExtension_bestScore = score;
            start = position1;
        }

        position1--; position2--;
    }

	// Correct for extra decrement
	start++;

    printf("\n%d to %d Dist=%d Score=%d\n", offset1, start, offset1 - start, ungappedExtension_bestScore);
/*
	// Starting at right/last hit position again
	queryPosition = queryHit + 1;
	subjectPosition = subjectEnd = subjectHit + 1;
	changeSinceBest = 0;

    // May need to alter dropoff so we also dropoff if below zero
    if (-ungappedExtension_bestScore > originalDropoff)
    {
    	dropoff = -ungappedExtension_bestScore;
    }

	// Extend end of alignment until dropoff
	while (changeSinceBest > dropoff)
	{
		changeSinceBest += (*queryPosition)[*subjectPosition];
		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				queryPosition++; subjectPosition++;
				changeSinceBest = (*queryPosition)[*subjectPosition];
			}
			while (changeSinceBest > 0);

			subjectEnd = subjectPosition;

			// Check need for change in dropoff
            if ((dropoff = -ungappedExtension_bestScore) < originalDropoff)
            {
            	dropoff = originalDropoff;
            }
        }
		queryPosition++; subjectPosition++;
	}

	// Correct for extra increment
	subjectEnd--;

    // Record the point we got to extending forwards
    ungappedExtension_subjectEndReached = subjectPosition;

    // If extension scored above trigger for gapping, create object and return it
    if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger)
    {
    	int4 diagonal;
        struct ungappedExtension* newUngappedExtension;
        newUngappedExtension = memBlocks_newEntry(ungappedExtension_extensions);

        // Calculate diagonal
        diagonal = (subjectHit - subject) - (queryHit - PSSMatrix.matrix);

        // Determine offsets from pointers
        newUngappedExtension->start.subjectOffset = subjectStart - subject;
        newUngappedExtension->end.subjectOffset = subjectEnd - subject;
        newUngappedExtension->start.queryOffset = newUngappedExtension->start.subjectOffset - diagonal;
        newUngappedExtension->end.queryOffset = newUngappedExtension->end.subjectOffset - diagonal;

        // Find the seed point
        newUngappedExtension->seed = ungappedExtension_findProteinSeed(newUngappedExtension, PSSMatrix, subject);
        // Initialize next to null
        newUngappedExtension->next = NULL;
        newUngappedExtension->nominalScore = ungappedExtension_bestScore;
        newUngappedExtension->status = ungappedExtension_UNGAPPED;

        return newUngappedExtension;
    }
    else
    {
    	return NULL;
    }*/
    return NULL;
}




// Out of date code
// ----------------

/*
// Given start and end query/subject positions, finds the highest scoring
// ungapped region between them and extending outwards until dropoff
// (Used to score CAFE frames)
inline void ungappedExtension_localExtend()
{
	int2** queryPosition, **newQueryStart = NULL, **newQueryEnd = NULL;
	unsigned char* subjectPosition, *newSubjectStart = NULL, *newSubjectEnd = NULL;
	int4 score = 0, bestScore = 0;
	int4 changeSinceBest = 0;
	int4 dropoff, originalDropoff;

    originalDropoff = dropoff = -statistics_ungappedNominalDropoff;

	// Start at queryEnd,subjectEnd (end of cafe frame)
	queryPosition = ungappedExtension_queryEnd;
	subjectPosition = ungappedExtension_subjectEnd;

	// Move through the frame backwards finding the highest scoring region
	while (queryPosition >= ungappedExtension_queryStart)
	{
    	score += (*queryPosition)[*subjectPosition];

        // Reset negative values to zero, recording new end position
        if (score < 0)
        {
        	ungappedExtension_queryEnd = queryPosition;
            ungappedExtension_subjectEnd = subjectPosition;
        	score = 0;
        }
        // Best score yet
        else if (score > bestScore)
        {
         	// Record start and end of high scoring region
			newQueryEnd = ungappedExtension_queryEnd;
			newSubjectEnd = ungappedExtension_subjectEnd;
			newQueryStart = queryPosition;
			newSubjectStart = subjectPosition;
            bestScore = score;
        }

        queryPosition--; subjectPosition--;
	}

    // If no region scored above zero, return zero
    if (newQueryStart == NULL)
	{
		ungappedExtension_nominalScore = 0;
		return;
	}

    // Else mark start and end of best region
    ungappedExtension_queryStart = newQueryStart;
	ungappedExtension_subjectStart = newSubjectStart;
    ungappedExtension_queryEnd = newQueryEnd;
	ungappedExtension_subjectEnd = newSubjectEnd;

	// Extend end forwards until dropoff
	queryPosition = ungappedExtension_queryEnd;
	subjectPosition = ungappedExtension_subjectEnd;
    changeSinceBest = 0;

	// Extend end of alignment until dropoff
	while (changeSinceBest > dropoff)
	{
		changeSinceBest += (*queryPosition)[*subjectPosition];
		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				bestScore += changeSinceBest;
				queryPosition++; subjectPosition++;
				changeSinceBest = (*queryPosition)[*subjectPosition];
			}
			while (changeSinceBest > 0);

			ungappedExtension_queryEnd = queryPosition;
			ungappedExtension_subjectEnd = subjectPosition;
        }
		queryPosition++; subjectPosition++;
	}

	// Correct for extra increment
	ungappedExtension_queryEnd--;
	ungappedExtension_subjectEnd--;

    // Extend the start of the hit forwards until dropoff
	queryPosition = ungappedExtension_queryStart - 1;
	subjectPosition = ungappedExtension_subjectStart - 1;
	while (changeSinceBest > dropoff)
	{
		changeSinceBest += (*queryPosition)[*subjectPosition];
		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				bestScore += changeSinceBest;
				queryPosition--; subjectPosition--;
				changeSinceBest = (*queryPosition)[*subjectPosition];
			}
			while (changeSinceBest > 0);

			ungappedExtension_queryStart = queryPosition;
			ungappedExtension_subjectStart = subjectPosition;
		}
		queryPosition--; subjectPosition--;
	}

	// Correct for extra decrement
	ungappedExtension_queryStart++;
	ungappedExtension_subjectStart++;

	ungappedExtension_nominalScore = bestScore;
}

void ungappedExtension_findBest(struct PSSMatrix PSSMatrix, unsigned char* subject, int4 subjectLength)
{
	int2 **queryPosition, **bestQueryPosition;
	unsigned char *subjectPosition, *bestSubjectPosition;
    int4 score, bestScore, offset, match;

    bestScore = 0;
    bestQueryPosition = PSSMatrix.matrix;
    bestSubjectPosition = subject;

    // For each diagonal starting with query = 0
    offset = 0;
    while (offset < subjectLength)
    {
    	// Start at beginning of diagonal
	    queryPosition = PSSMatrix.matrix;
		subjectPosition = subject + offset;

        // Move along until end of diagonal (hit sentinal byte)
        score = 0;
        while ((match = (*queryPosition)[*subjectPosition]) > constants_sentinalScore)
		{
        	// Accumulate scores
			score += match;
            // If score drops below zero, reset it
            if (score < 0)
            {
            	score = 0;
            }
            // Check for best score yet;
            else if (score > bestScore)
            {
            	bestQueryPosition = queryPosition;
                bestSubjectPosition = subjectPosition;
                bestScore = score;
            }

            queryPosition++; subjectPosition++;
        }

    	offset++;
    }

    // For each diagonal starting with subject = 0
    offset = 0;
    while (offset < PSSMatrix.length)
    {
    	// Start at beginning of diagonal
	    queryPosition = PSSMatrix.matrix + offset;
		subjectPosition = subject;

        // Move along until end of diagonal (hit sentinal byte)
        score = 0;
        while ((match = (*queryPosition)[*subjectPosition]) > constants_sentinalScore)
		{
        	// Accumulate scores
			score += match;
            // If score drops below zero, reset it
            if (score < 0)
            {
            	score = 0;
            }
            // Check for best score yet;
            else if (score > bestScore)
            {
            	bestQueryPosition = queryPosition;
                bestSubjectPosition = subjectPosition;
                bestScore = score;
            }

            queryPosition++; subjectPosition++;
        }

    	offset++;
    }

//    printf("Best=%d,%d Score=%d\n", bestQueryPosition - PSSMatrix.matrix,
//                                    bestSubjectPosition - subject, bestScore);
}

// Returns a copy of the current ungappedExtension information bundled into a structure
struct ungappedExtension* ungappedExtension_copy(struct PSSMatrix PSSMatrix, unsigned char* subject)
{
	struct ungappedExtension* newUngappedExtension;

	newUngappedExtension = memBlocks_newEntry(ungappedExtension_extensions);

	// Copy all of the query/subject start/end pointer information
	newUngappedExtension->start.queryPosition = ungappedExtension_queryStart;
	newUngappedExtension->end.queryPosition = ungappedExtension_queryEnd;
	newUngappedExtension->start.subjectPosition = ungappedExtension_subjectStart;
	newUngappedExtension->end.subjectPosition = ungappedExtension_subjectEnd;
	// Determine offsets from pointers
	newUngappedExtension->start.queryOffset = ungappedExtension_queryStart - PSSMatrix.matrix;
	newUngappedExtension->end.queryOffset = ungappedExtension_queryEnd - PSSMatrix.matrix;
	newUngappedExtension->start.subjectOffset = ungappedExtension_subjectStart - subject;
	newUngappedExtension->end.subjectOffset = ungappedExtension_subjectEnd - subject;

	newUngappedExtension->nominalScore = ungappedExtension_nominalScore;
	newUngappedExtension->scoringStatus = ungappedExtension_UNGAPPED;

	// Find the seed point
	newUngappedExtension->seed = ungappedExtension_findSeed(newUngappedExtension, PSSMatrix, subject);

	// Initialize next (for building a linked list of extensions) to null
	newUngappedExtension->next = NULL;

	return newUngappedExtension;
}

// Copy the list of "from" ungappedExtensions to "to" and clear their scores and status in the
// process
struct ungappedExtension* ungappedExtension_duplicate(struct ungappedExtension* from,
                          struct PSSMatrix PSSMatrix, unsigned char* subject)
{
	struct ungappedExtension *newUngappedExtension, *newList = NULL;

    while (from != NULL)
    {
    	// Create new extension
        newUngappedExtension = memBlocks_newEntry(ungappedExtension_extensions);

//        printf("Duplicate: %d,%d -> %d,%d\n", from->start.queryOffset, from->start.subjectOffset,
//               from->end.queryOffset, from->end.subjectOffset); fflush(stdout);

        // Copy details
        newUngappedExtension->start.queryOffset = from->start.queryOffset;
        newUngappedExtension->end.queryOffset = from->end.queryOffset;
        newUngappedExtension->start.subjectOffset = from->start.subjectOffset;
        newUngappedExtension->end.subjectOffset = from->end.subjectOffset;
        // Determine pointers from offset
        newUngappedExtension->start.queryPosition
        	= newUngappedExtension->start.queryOffset + PSSMatrix.matrix;
        newUngappedExtension->end.queryPosition
        	= newUngappedExtension->end.queryOffset + PSSMatrix.matrix;
        newUngappedExtension->start.subjectPosition
        	= newUngappedExtension->start.subjectOffset + subject;
        newUngappedExtension->end.subjectPosition
        	= newUngappedExtension->end.subjectOffset + subject;

//        newUngappedExtension->offset = from->offset;
        newUngappedExtension->nominalScore = 0;
        newUngappedExtension->scoringStatus = ungappedExtension_UNGAPPED;

        newUngappedExtension->seed.queryOffset = from->seed.queryOffset;
        newUngappedExtension->seed.subjectOffset = from->seed.subjectOffset;
        newUngappedExtension->seed.queryPosition = from->seed.queryOffset + PSSMatrix.matrix;
        newUngappedExtension->seed.subjectPosition = from->seed.subjectOffset + subject;

        // Add to the front of the list
		newUngappedExtension->next = newList;
        newList = newUngappedExtension;

        // Advance to next extension
        from = from->next;
    }

    return newList;
}

uint4 readLittleEndianInt(unsigned char* address)
{
	uint4 value = 0;
    unsigned char* a = address + 3;

    while (a >= address)
    {
    	value *= 256;
		value += *a;
        a--;
    }

    return value;
}

// Returns a copy of an ungappedExtension constructed using given hspindex information
struct ungappedExtension* ungappedExtension_construct(struct PSSMatrix PSSMatrix,
                                     unsigned char* subject, int4* hspInformation)
{
	struct ungappedExtension* newUngappedExtension;

	newUngappedExtension = memBlocks_newEntry(ungappedExtension_extensions);

	// Copy all of the query/subject start/end pointer information
	newUngappedExtension->start.queryOffset = readLittleEndianInt((unsigned char*)hspInformation); hspInformation++;
	newUngappedExtension->start.subjectOffset = readLittleEndianInt((unsigned char*)hspInformation); hspInformation++;
	newUngappedExtension->end.queryOffset = readLittleEndianInt((unsigned char*)hspInformation); hspInformation++;
	newUngappedExtension->end.subjectOffset = readLittleEndianInt((unsigned char*)hspInformation); hspInformation++;

    // Determine offsets from pointers
	newUngappedExtension->start.queryPosition
    	= newUngappedExtension->start.queryOffset + PSSMatrix.matrix;
	newUngappedExtension->end.queryPosition
    	= newUngappedExtension->end.queryOffset + PSSMatrix.matrix;
	newUngappedExtension->start.subjectPosition
    	= newUngappedExtension->start.subjectOffset + subject;
	newUngappedExtension->end.subjectPosition
    	= newUngappedExtension->end.subjectOffset + subject;

//	newUngappedExtension->offset
//    	= newUngappedExtension->end.subjectOffset - newUngappedExtension->end.queryOffset;
	newUngappedExtension->nominalScore = readLittleEndianInt((unsigned char*)hspInformation); hspInformation++;
	newUngappedExtension->scoringStatus = ungappedExtension_UNGAPPED;

	// Find the seed point
	newUngappedExtension->seed = ungappedExtension_findSeed(newUngappedExtension, PSSMatrix, subject);

	// Initialize next (for building a linked list of extensions) to null
	newUngappedExtension->next = NULL;

	return newUngappedExtension;
}
*/
