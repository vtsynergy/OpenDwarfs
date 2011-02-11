// search.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Functions for performing stages 1 & 2 of the search algorithm

#include "blast.h"
#include "wordLookupDFA.h"
extern unsigned char * wordLookupDFA;
extern struct groupFP *wordLookupDFA_groupsFP;

// Shucai
// Search a protein database using 1-hit extension mode
void search_protein1hit(struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, struct sequenceData* sequenceData,
						uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
	unsigned char *subject, *sequenceEnd, *address;
	int4 subjectLength, subjectOffset, wordLengthMinusOne, count = 0;
	unsigned char currentWord, *currentBlock;
    struct group* currentGroup;
	//Shucai
	struct groupFP *currentGroupFP;
	unsigned char *startAddressFP = readdb_sequences;

	uint2* queryOffsets, queryOffset;
	struct ungappedExtension* ungappedExtension;
    int4 diagonal;

	//Shucai
    unsigned char** lastHit;
	uint4 *lastHitFP;

    wordLengthMinusOne = parameters_wordSize - 1;

    while (sequenceCount < numSequences)
    {
    	descriptionLength = sequenceData[sequenceCount].descriptionLength;
    	descriptionStart = sequenceData[sequenceCount].descriptionStart;
    	subjectLength = sequenceData[sequenceCount].sequenceLength;
    	encodedLength = sequenceData[sequenceCount].encodedLength;
		address = subject = sequenceData[sequenceCount].sequence;
		
        // New sequence, new possible alignment
        alignments_currentAlignment = NULL;

        // Only process sequence if at least as long as the word length
        if (subjectLength >= parameters_wordSize)
        {
            // Start at 000 state in Finite State Automata
			//Shucai
            //currentGroup = wordLookupDFA_groups;
			currentGroupFP = wordLookupDFA_groupsFP;

            // Read first wordLength - 1 chars and advance
            count = wordLengthMinusOne;
            while (count > 0)
            {
                if (*address < wordLookupDFA_numCodes)
				{
					//Shucai
                    //currentGroup = currentGroup->nextGroups + *address;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
                }
				else
				{
					//Shucai
                    //currentGroup = currentGroup->nextGroups;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
				}
                address++;
                count--;
            }

            // Read the rest of the codes, using the Finite State Automata defined
            // by wordLookupDFA to get the query positions of int4erest
            sequenceEnd = subject + subjectLength;
            while (address < sequenceEnd)
            {
				//Shucai
                //currentBlock = currentGroup->nextWords;
				currentBlock = &wordLookupDFA[currentGroupFP->nextWords];

                // If current code is a regular letter
                if (*address < wordLookupDFA_numCodes)
                {
                    // Use it
                    currentWord = currentBlock[*address];
					//Shucai
                    //currentGroup = currentGroup->nextGroups + *address;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
                }
                else
                {
                    // Else check if we've reached end of the file
                    if (address >= sequenceEnd)
                        break;

                    // If not, we've read a wild code. Use first code instead
                    currentWord = *currentBlock;
					//Shucai
                    //currentGroup = currentGroup->nextGroups;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
                }

                if (currentWord)
                {
                    // Calculate subject offset
                    subjectOffset = address - subject;

                    // At least one query position, stored at an extenal address
                    queryOffsets = ((uint2*)currentBlock) - currentWord;

                    // If the zero flag is stored at the first query position
                    if (!*queryOffsets)
                    {
                        // Go to an outside address for additional positions
                        queryOffsets = wordLookupDFA_additionalQueryPositions
                                     + (*(queryOffsets + 1) * constants_max_int2) + *(queryOffsets + 2);
                    }

                    do
                    {
                        queryOffset = *queryOffsets;

                        #ifndef NO_STAGE2
                        // Calculate the diagonal this hit is on
                        diagonal = subjectOffset - queryOffset;

                        // If we have not extended past this point on this diagonal
                        //Shucai
						//lastHit = hitMatrix_furthest + diagonal;
						lastHitFP = hitMatrix_furthestFP + diagonal;
						//Shucai
						//if (*lastHit < address)
                        if (*lastHitFP < address - startAddressFP)
                        {
                            // Increment tally number of extensions
                            blast_numUngappedExtensions++;

                            // If only one hit triggered this extension
							// Shucai
                            ungappedExtension
                            	//= ungappedExtension_oneHitExtend(PSSMatrix.matrix + queryOffset,
								= ungappedExtension_oneHitExtend(queryOffset,
                                    address, PSSMatrix, PSSMatrixFP, subject, startAddressFP);

                            // Update furthest reached value for the diagonal
							// Shucai
                            //*lastHit = ungappedExtension_subjectEndReached;
							*lastHitFP = ungappedExtension_subjectEndReachedFP;
                            // If extension scores high enough to trigger gapping
                            if (ungappedExtension)
                            {
                                // Increment count of number of trigger extensions
                                blast_numTriggerExtensions++;

                                // Create new alignment if needed
                                if (alignments_currentAlignment == NULL)
                                {
                                    // Create new alignment object
                                    alignments_createNew(descriptionStart, descriptionLength,
                                                         subject, subjectLength, encodedLength);
                                }

                                // Add this extension to the alignment
                                alignments_addUngappedExtension(ungappedExtension);
                            }
                        }
                        #endif

                        queryOffsets++; blast_numHits++;
                    }
                    while (*queryOffsets);
                }

                address++;
            }
        }

        sequenceCount++;

        // Every so often print status .
        if ((sequenceCount % tickFrequency == 0) &&
            parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
        {
            #ifndef VERBOSE
            printf(".");
            fflush(stdout);
			#endif
        }
    }
}

// Search a protein database using 2-hit extension mode
void search_protein2hit(struct PSSMatrix PSSMatrix, struct PSSMatrixFP PSSMatrixFP, struct sequenceData* sequenceData,
						uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
	unsigned char *subject, *sequenceEnd;
	int4 subjectLength, subjectOffset, wordLengthMinusOne, count = 0;
	unsigned char currentWord, *currentBlock, *address;
    struct group* currentGroup;
	//Shucai
	struct groupFP *currentGroupFP;
	unsigned char *startAddressFP = readdb_sequences;

	uint2* queryOffsets, queryOffset;
	struct ungappedExtension* ungappedExtension;
    int4 diagonal, distance;
	//Shucai
    //unsigned char** lastHit;
	uint4 *lastHitFP;

    wordLengthMinusOne = parameters_wordSize - 1;

    while (sequenceCount < numSequences)
    {
    	descriptionLength = sequenceData[sequenceCount].descriptionLength;
    	descriptionStart = sequenceData[sequenceCount].descriptionStart;
    	subjectLength = sequenceData[sequenceCount].sequenceLength;
    	encodedLength = sequenceData[sequenceCount].encodedLength;
		address = subject = sequenceData[sequenceCount].sequence;

        // New sequence, new possible alignment
        alignments_currentAlignment = NULL;

        // Only process sequence if at least as long as the word length
        if (subjectLength >= parameters_wordSize)
        {
            // Start at 000 state in Finite State Automata
			// Shucai
            //currentGroup = wordLookupDFA_groups;
			currentGroupFP = wordLookupDFA_groupsFP;

            // Read first wordLength - 1 chars and advance
            count = wordLengthMinusOne;
            while (count > 0)
            {
                if (*address < wordLookupDFA_numCodes)
				{
					//Shucai
                    //currentGroup = currentGroup->nextGroups + *address;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
                }
				else
				{
					//Shucai
                    //currentGroup = currentGroup->nextGroups;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
				}

                address++;
                count--;
            }

            // Read the rest of the codes, using the Finite State Automata defined
            // by wordLookupDFA to get the query positions of int4erest
            sequenceEnd = subject + subjectLength;
            while (address < sequenceEnd)
            {
				//Shucai
                //currentBlock = currentGroup->nextWords;
				currentBlock = &wordLookupDFA[currentGroupFP->nextWords];

                // If current code is a regular letter
                if (*address < wordLookupDFA_numCodes)
                {
                    // Use it
                    currentWord = currentBlock[*address];
					//Shucai
                    //currentGroup = currentGroup->nextGroups + *address;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
                }
                else
                {
                    // Else check if we've reached end of the file
                    if (address >= sequenceEnd)
                        break;

                    // If not, we've read a wild code. Use first code instead
                    currentWord = currentBlock[0];
					//Shucai
                    //currentGroup = currentGroup->nextGroups;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
                }

                if (currentWord)
                {
                    // Calculate subject offset
                    subjectOffset = address - subject;

                    // If at least one query position, stored at an extenal address
                    queryOffsets = ((uint2*)currentBlock) - currentWord;

                    // If the zero flag is stored at the first query position
                    if (!*queryOffsets)
                    {
                        // Go to an outside address for additional positions
                        queryOffsets = wordLookupDFA_additionalQueryPositions
                                     + (*(queryOffsets + 1) * constants_max_int2) + *(queryOffsets + 2);
                    }

                    do
                    {
                        queryOffset = *queryOffsets;

                        #ifndef NO_STAGE2
                        // Calculate the diagonal this hit is on
                        diagonal = subjectOffset - queryOffset;

                        // Calculate distance since last hit
                        //Shucai
						//lastHit = hitMatrix_furthest + diagonal;
						lastHitFP = hitMatrix_furthestFP + diagonal;
                        //Shucai
						//distance = address - *lastHit;
						distance = (address - startAddressFP) - *lastHitFP;

						//Shucai


                        if (distance >= parameters_A)
                        {
                            // Too far apart, update furthest
							//Shucai
                            //*lastHit = address;
							*lastHitFP = address - startAddressFP;
                        }
                        else if (distance >= parameters_overlap)
                        {
                            // Not overlaping - extension triggered
                            // Increment tally number of extensions
                            blast_numUngappedExtensions++;

                            // Perform ungapped extension start between query/subject start/end
                            // and extending outwards in each direction
                            ungappedExtension
								//= ungappedExtension_extend(PSSMatrix.matrix + queryOffset,
                            	= ungappedExtension_extend(queryOffset,
								//Shucai
                                     //address, *lastHit, PSSMatrix, PSSMatrixFP, subject);
                                     address, *lastHitFP, PSSMatrix, PSSMatrixFP, subject, startAddressFP);

                            // Update furthest reached value for the diagonal
                            //Shucai
							//*lastHit = ungappedExtension_subjectEndReached;
							*lastHitFP = ungappedExtension_subjectEndReachedFP;

                            // If extension scores high enough to trigger gapping
                            if (ungappedExtension)
                            {
                                // Increment count of number of trigger extensions
                                blast_numTriggerExtensions++;

                                // Create new alignment if needed
                                if (alignments_currentAlignment == NULL)
                                {
                                    // Create new alignment object using subject with wilds
                                    alignments_createNew(descriptionStart, descriptionLength, subject,
                                                         subjectLength, encodedLength);
                                }

                                // Add this extension to the alignment
                                alignments_addUngappedExtension(ungappedExtension);
                            }
                        }
                        #endif

                        queryOffsets++; blast_numHits++;
                    }
                    while (*queryOffsets);
                }

                address++;
            }
        }

        sequenceCount++;

        // Every so often print status .
        if ((sequenceCount % tickFrequency == 0) &&
            parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
        {
            #ifndef VERBOSE
            printf(".");
            fflush(stdout);
			#endif
        }
    }
}

// Search a nucleotide database using 1-hit extension mode
void search_nucleotide(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                       uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0, codeword;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    uint4 numPackedBytes, numRemaining;
	unsigned char *subject, *sequenceEnd, prevByte, nextByte, rightMatches, leftMatches;
    unsigned char previousCode, *subjectPosition, *address;
	int4 subjectLength, subjectOffset, wordLengthMinusOne;
	struct ungappedExtension* ungappedExtension;
	int2 *queryOffsets, tableEntry, queryOffset, singleQueryOffset[2];
    int4 previousByteDistance;
    int4 diagonal;
    unsigned char** lastHit;

    wordLengthMinusOne = parameters_wordSize - 1;
	previousByteDistance = (parameters_wordTableLetters + 4);
	singleQueryOffset[1] = 0;

    while (sequenceCount < numSequences)
    {
    	descriptionLength = sequenceData[sequenceCount].descriptionLength;
    	descriptionStart = sequenceData[sequenceCount].descriptionStart;
    	subjectLength = sequenceData[sequenceCount].sequenceLength;
    	encodedLength = sequenceData[sequenceCount].encodedLength;
		address = subject = sequenceData[sequenceCount].sequence;

        // New sequence, new possible alignment
        alignments_currentAlignment = NULL;

        // Determine number of packed bytes and number of remaining letters
        numPackedBytes = subjectLength / 4;
		numRemaining = subjectLength % 4;

        // Only process sequence if at least as long as the word length
        if (subjectLength >= parameters_wordSize)
        {
            sequenceEnd = subject + numPackedBytes;

            // Read first char and advance
            previousCode = *address;
            address++;

            // Traverse until end of sequence
            while (address < sequenceEnd)
            {
            	// Read next char and update codeword
            	codeword = *address | (previousCode << 8);
                previousCode = *address;
				tableEntry = nucleotideLookup[codeword];

                // Calculate subject offset
                subjectOffset = address - subject;

                #ifdef BITLOOKUP
				if (nucleotideLookup_bitLookup[codeword >> 5] &
                    (1 << (codeword & nucleotideLookup_mask)))
				#endif
                if (tableEntry)
                {
                	if (tableEntry > 0)
                    {
                    	// Just one query position
                    	singleQueryOffset[0] = queryOffset = tableEntry;
						queryOffsets = singleQueryOffset;
                    }
                    else
                    {
                        // Multiple query positions
						queryOffsets = nucleotideLookup_additionalPositions - tableEntry;
                        queryOffset = *queryOffsets;
                    }

                    subjectPosition = address + 1;
                    do
                    {
                        // Calculate number of matches to left and right of hit
                        nextByte = PSSMatrix.bytePackedCodes[queryOffset];
                        rightMatches = PSSMatrix_packedRightMatches[nextByte ^ *subjectPosition];
                        prevByte = PSSMatrix.bytePackedCodes[queryOffset - previousByteDistance];
                        leftMatches = PSSMatrix_packedLeftMatches[prevByte ^ *(address - parameters_wordTableBytes)];

                        if (rightMatches + leftMatches >= parameters_wordExtraLetters)
                        {
                            #ifndef NO_STAGE2
                            // Calculate subject offset
                            subjectOffset = subjectPosition - subject;

                            // Calculate the diagonal this hit is on
                            diagonal = (subjectOffset * 4 - queryOffset + hitMatrix_queryLength)
                                     & hitMatrix_diagonalMask;

                            // If we have not extended past this point on this diagonal
                            lastHit = hitMatrix_furthest + diagonal;

                            #ifdef VERBOSE
                            if (parameters_verboseDloc == descriptionStart)
                                printf("Hit %d,%d\n", queryOffset, subjectOffset * 4);
                            #endif

                            if (*lastHit < address)
                            {
                                // Perform ungapped extension
                                ungappedExtension
                                	= ungappedExtension_nucleotideExtend(queryOffset,
                                      subjectOffset, PSSMatrix, subject, subjectLength);

                                // Update furthest reached value for the diagonal
                                *lastHit = ungappedExtension_subjectEndReached;

                                #ifdef VERBOSE
                                if (parameters_verboseDloc == descriptionStart)
                                    printf("UngappedExtension %d,%d Score=%d\n", queryOffset,
                                    subjectOffset * 4, ungappedExtension_bestScore);
                                if (parameters_verboseDloc == descriptionStart && ungappedExtension)
                                	ungappedExtension_print(ungappedExtension);
                                #endif

                                // If extension scores high enough to trigger gapping
                                if (ungappedExtension)
                                {
                                    // Increment count of number of trigger extensions
                                    blast_numTriggerExtensions++;

                                    // Create new alignment if needed
                                    if (alignments_currentAlignment == NULL)
                                    {
                                        alignments_createNew(descriptionStart, descriptionLength, subject,
                                                             subjectLength, encodedLength);
                                    }

                                    // Add this extension to the alignment
                                    alignments_addUngappedExtension(ungappedExtension);
                                }

                                blast_numUngappedExtensions++;
                            }
                            #endif
                            blast_numHits++;
                        }

                        queryOffsets++;
                        queryOffset = *queryOffsets;
					}
                    while (queryOffset);
                }

                address++;
            }
        }

        sequenceCount++;

        // Every so often print status .
        if ((sequenceCount % tickFrequency == 0) &&
            parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
        {
            #ifndef VERBOSE
            printf(".");
            fflush(stdout);
            #endif
        }
    }
}

// Search a nucleotide database using 1-hit extension mode with large wordsize > 14
void search_nucleotide_longWord(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                                uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0, codeword;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    uint4 numPackedBytes, numRemaining;
	unsigned char *subject, *sequenceEnd, prevByte, nextByte, rightMatches, leftMatches;
    unsigned char previousCode, *subjectPosition, *address;
	int4 subjectLength, subjectOffset, wordLengthMinusOne;
	struct ungappedExtension* ungappedExtension;
	int2 *queryOffsets, tableEntry, queryOffset, singleQueryOffset[2];
    int4 previousByteDistance;
    int4 diagonal, extraBytesNeeded;
    unsigned char** lastHit;

    wordLengthMinusOne = parameters_wordSize - 1;
	previousByteDistance = parameters_wordTableLetters + parameters_wordExtraBytes * 4 + 4;
	singleQueryOffset[1] = 0;

    while (sequenceCount < numSequences)
    {
    	descriptionLength = sequenceData[sequenceCount].descriptionLength;
    	descriptionStart = sequenceData[sequenceCount].descriptionStart;
    	subjectLength = sequenceData[sequenceCount].sequenceLength;
    	encodedLength = sequenceData[sequenceCount].encodedLength;
		address = subject = sequenceData[sequenceCount].sequence;

        // Determine number of packed bytes and number of remaining letters
        numPackedBytes = subjectLength / 4;
		numRemaining = subjectLength % 4;

        // Only process sequence if at least as long as the word length
        if (subjectLength >= parameters_wordSize)
        {
            sequenceEnd = subject + numPackedBytes;

            // Read first char and advance
            previousCode = *address;
            address++;

            // Traverse until end of sequence
            while (address < sequenceEnd)
            {
            	// Read next char and update codeword
            	codeword = *address | (previousCode << 8);
                previousCode = *address;
				tableEntry = nucleotideLookup[codeword];

                // Calculate subject offset
                subjectOffset = address - subject;

                #ifdef BITLOOKUP
				if (nucleotideLookup_bitLookup[codeword >> 5] &
                    (1 << (codeword & nucleotideLookup_mask)))
				#endif
                if (tableEntry)
                {
                	if (tableEntry > 0)
                    {
                    	// Just one query position
                    	singleQueryOffset[0] = queryOffset = tableEntry;
						queryOffsets = singleQueryOffset;
                    }
                    else
                    {
                        // Multiple query positions
						queryOffsets = nucleotideLookup_additionalPositions - tableEntry;
                        queryOffset = *queryOffsets;
                    }

                    subjectPosition = address + 1;
                    do
                    {
						extraBytesNeeded = parameters_wordExtraBytes;

						while (extraBytesNeeded)
						{
							// Check for matching bytes to right
                            if (*subjectPosition != PSSMatrix.bytePackedCodes[queryOffset])
                            	break;

//							printf("Match %d,%d\n", queryOffset, (subjectPosition - subject)*4);

                            extraBytesNeeded--;
                            subjectPosition++;
                            subjectOffset++;
                            queryOffset+=4;
						}

						if (!extraBytesNeeded)
						{
                            // Calculate number of matches to left and right of hit
                            nextByte = PSSMatrix.bytePackedCodes[queryOffset];
                            rightMatches = PSSMatrix_packedRightMatches[nextByte ^ *(subjectPosition)];
                            prevByte = PSSMatrix.bytePackedCodes[queryOffset - previousByteDistance];
                            leftMatches = PSSMatrix_packedLeftMatches[prevByte ^ *(address - parameters_wordTableBytes)];

//                            printf("prev at %d,%d\n", queryOffset - previousByteDistance);
//                        	printf("part matches=[%d,%d]\n", leftMatches, rightMatches);

                            if (rightMatches + leftMatches >= parameters_wordExtraLetters)
                            {
//                            	printf("Hit! dloc=%d\n", descriptionStart);
                                #ifndef NO_STAGE2
                                // Calculate subject offset
                                subjectOffset = subjectPosition - subject;

                                // Calculate the diagonal this hit is on
                                diagonal = (subjectOffset * 4 - queryOffset + hitMatrix_queryLength)
                                        & hitMatrix_diagonalMask;

                                // If we have not extended past this point on this diagonal
                                lastHit = hitMatrix_furthest + diagonal;

                                #ifdef VERBOSE
                                if (parameters_verboseDloc == descriptionStart)
                                    printf("Hit %d,%d\n", queryOffset, subjectOffset * 4);
                                #endif
                                if (*lastHit < address)
                                {
                                    // Perform ungapped extension
                                    ungappedExtension
                                        = ungappedExtension_nucleotideExtend(queryOffset,
                                        subjectOffset, PSSMatrix, subject, subjectLength);

                                    // Update furthest reached value for the diagonal
                                    *lastHit = ungappedExtension_subjectEndReached;

                                    #ifdef VERBOSE
                                    if (parameters_verboseDloc == descriptionStart)
                                        printf("UngappedExtension %d,%d Score=%d\n", queryOffset,
                                        subjectOffset * 4, ungappedExtension_bestScore);
                                    if (parameters_verboseDloc == descriptionStart && ungappedExtension)
                                        ungappedExtension_print(ungappedExtension);
                                    #endif

                                    // If extension scores high enough to trigger gapping
                                    if (ungappedExtension)
                                    {
                                        // Increment count of number of trigger extensions
                                        blast_numTriggerExtensions++;

                                        // Create new alignment if needed
                                        if (alignments_currentAlignment == NULL)
                                        {
                                            alignments_createNew(descriptionStart, descriptionLength, subject,
                                                                subjectLength, encodedLength);
                                        }

                                        // Add this extension to the alignment
                                        alignments_addUngappedExtension(ungappedExtension);
                                    }

                                    blast_numUngappedExtensions++;
                                }
                                #endif
                                blast_numHits++;
                            }
                        }

                        queryOffsets++;
                        queryOffset = *queryOffsets;
					}
                    while (queryOffset);
                }

                address++;
            }
        }

        sequenceCount++;

        // Every so often print status .
        if ((sequenceCount % tickFrequency == 0) &&
            parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
        {
            #ifndef VERBOSE
            printf(".");
            fflush(stdout);
            #endif
        }

        // New sequence, new possible alignment
        alignments_currentAlignment = NULL;
    }
}

// Search a nucleotide database using 1-hit extension mode, using a large word lookup table
// due to long query sequence
void search_nucleotide_largeTable(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                                  uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0, codeword;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    uint4 numPackedBytes, numRemaining;
	unsigned char *subject, *sequenceEnd, prevByte, nextByte, rightMatches, leftMatches;
    unsigned char previousCode, *subjectPosition, *address;
	int4 subjectLength, subjectOffset, wordLengthMinusOne;
	struct ungappedExtension* ungappedExtension;
	int4 *queryOffsets, tableEntry, queryOffset, singleQueryOffset[2];
    int4 previousByteDistance;
    int4 diagonal;
    unsigned char** lastHit;

    wordLengthMinusOne = parameters_wordSize - 1;
	previousByteDistance = (parameters_wordTableLetters + 4);
	singleQueryOffset[1] = 0;

    while (sequenceCount < numSequences)
    {
    	descriptionLength = sequenceData[sequenceCount].descriptionLength;
    	descriptionStart = sequenceData[sequenceCount].descriptionStart;
    	subjectLength = sequenceData[sequenceCount].sequenceLength;
    	encodedLength = sequenceData[sequenceCount].encodedLength;
		address = subject = sequenceData[sequenceCount].sequence;

        // New sequence, new possible alignment
        alignments_currentAlignment = NULL;

        // Determine number of packed bytes and number of remaining letters
        numPackedBytes = subjectLength / 4;
		numRemaining = subjectLength % 4;

        // Only process sequence if at least as long as the word length
        if (subjectLength >= parameters_wordSize)
        {
            sequenceEnd = subject + numPackedBytes;

            // Read first char and advance
            previousCode = *address;
            address++;

            // Traverse until end of sequence
            while (address < sequenceEnd)
            {
            	// Read next char and update codeword
            	codeword = *address | (previousCode << 8);
                previousCode = *address;
				tableEntry = nucleotideLookup_large[codeword];

                // Calculate subject offset
                subjectOffset = address - subject;

                #ifdef BITLOOKUP
				if (nucleotideLookup_bitLookup[codeword >> 5] &
                    (1 << (codeword & nucleotideLookup_mask)))
				#endif
                if (tableEntry)
                {
                	if (tableEntry > 0)
                    {
                    	// Just one query position
                    	singleQueryOffset[0] = queryOffset = tableEntry;
						queryOffsets = singleQueryOffset;
                    }
                    else
                    {
                        // Multiple query positions
						queryOffsets = nucleotideLookup_additionalPositions_large - tableEntry;
                        queryOffset = *queryOffsets;
                    }

                    subjectPosition = address + 1;
                    do
                    {
                        // Calculate number of matches to left and right of hit
                        nextByte = PSSMatrix.bytePackedCodes[queryOffset];
                        rightMatches = PSSMatrix_packedRightMatches[nextByte ^ *subjectPosition];
                        prevByte = PSSMatrix.bytePackedCodes[queryOffset - previousByteDistance];
                        leftMatches = PSSMatrix_packedLeftMatches[prevByte ^ *(address - parameters_wordTableBytes)];

                        if (rightMatches + leftMatches >= parameters_wordExtraLetters)
                        {
                            #ifndef NO_STAGE2
                            // Calculate subject offset
                            subjectOffset = subjectPosition - subject;

                            // Calculate the diagonal this hit is on
                            diagonal = (subjectOffset * 4 - queryOffset + hitMatrix_queryLength)
                                     & hitMatrix_diagonalMask;

                            // If we have not extended past this point on this diagonal
                            lastHit = hitMatrix_furthest + diagonal;

                            #ifdef VERBOSE
                            if (parameters_verboseDloc == descriptionStart)
                                printf("Hit %d,%d\n", queryOffset, subjectOffset * 4);
                            #endif
                            if (*lastHit < address)
                            {
                                // Perform ungapped extension
                                ungappedExtension
                                	= ungappedExtension_nucleotideExtend(queryOffset,
                                      subjectOffset, PSSMatrix, subject, subjectLength);

                                // Update furthest reached value for the diagonal
                                *lastHit = ungappedExtension_subjectEndReached;

                                #ifdef VERBOSE
                                if (parameters_verboseDloc == descriptionStart)
                                    printf("UngappedExtension %d,%d Score=%d\n", queryOffset,
                                    subjectOffset * 4, ungappedExtension_bestScore);
                                if (parameters_verboseDloc == descriptionStart && ungappedExtension)
                                	ungappedExtension_print(ungappedExtension);
                                #endif

                                // If extension scores high enough to trigger gapping
                                if (ungappedExtension)
                                {
                                    // Increment count of number of trigger extensions
                                    blast_numTriggerExtensions++;

                                    // Create new alignment if needed
                                    if (alignments_currentAlignment == NULL)
                                    {
                                        alignments_createNew(descriptionStart, descriptionLength, subject,
                                                             subjectLength, encodedLength);
                                    }

                                    // Add this extension to the alignment
                                    alignments_addUngappedExtension(ungappedExtension);
                                }

                                blast_numUngappedExtensions++;
                            }
                            #endif
                            blast_numHits++;
                        }

                        queryOffsets++;
                        queryOffset = *queryOffsets;
					}
                    while (queryOffset);
                }

                address++;
            }
        }

        sequenceCount++;

        // Every so often print status .
        if ((sequenceCount % tickFrequency == 0) &&
            parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
        {
            #ifndef VERBOSE
            printf(".");
            fflush(stdout);
            #endif
        }
    }
}

/*
// Search a nucleotide collection using an inverted index to speed-up BLAST search
void search_nucleotideIndex(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                            uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0, codeword;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
    uint4 numPackedBytes, numRemaining;
	unsigned char *subject, *sequenceEnd, prevByte, nextByte, rightMatches, leftMatches;
    unsigned char previousCode, *unpackedSubject;
	int4 subjectLength, subjectPosition, wordLengthMinusOne, lastSubjectPosition;
	struct ungappedExtension* ungappedExtension;
	int4 *queryPositions, tableEntry, queryPosition, singleQueryPosition[2], lastQueryPosition;
    int4 previousByteDistance;
    int4 diagonal;
	uint4 *sequencePositions, *descriptionLocations;
    unsigned char** lastHit;
    unsigned char* startAddress;
	struct indexCoordinate *coordinate;
    uint4 databaseOffset;
	uint4 time;
	unsigned char* indexFilename, *startIndexFile;
    struct readFile indexFile;
	int4 lastSequenceNumber;

    wordLengthMinusOne = parameters_wordSize - 1;
	previousByteDistance = (parameters_wordTableBytes * 4 + 4);
	singleQueryPosition[1] = 0;

	indexFilename = (char*)global_malloc(strlen(blast_searchDbFile) + 9);
	sprintf(indexFilename, "%s.index", blast_searchDbFile);

	// Open index file for reading, mapping contents to address
	indexFile = readFile_open(indexFilename);
	startIndexFile = (unsigned char*)indexFile.address;

    time = clock();
	startAddress = address;

    // Process the query and inverted index to find hits
    index_processQuery(startIndexFile, PSSMatrix, readdb_numberOfSequences);

    // TODO: Use size of volume instead of size of database

	sequencePositions = (uint4*)global_malloc(sizeof(uint4) * readdb_numberOfSequences);
	descriptionLocations = (uint4*)global_malloc(sizeof(uint4) * readdb_numberOfSequences);

    // Load description positions from disk
    lastSequenceNumber = -1;
    coordinate = index_getFirstCoordinate();
    while (coordinate != NULL)
    {
		if (coordinate->subjectNumber != lastSequenceNumber)
        {
        	sequencePositions[coordinate->subjectNumber]
            	= index_sequencePositions[coordinate->subjectNumber];
		}
        lastSequenceNumber = coordinate->subjectNumber;

        // Get next coordinate;
        coordinate = index_getNextCoordinate();
    }

    // Load sequence positions from disk
    lastSequenceNumber = -1;
    coordinate = index_getFirstCoordinate();
    while (coordinate != NULL)
    {
		if (coordinate->subjectNumber != lastSequenceNumber)
        {
			descriptionLocations[coordinate->subjectNumber]
            	= index_descriptionLocations[coordinate->subjectNumber];
		}
        lastSequenceNumber = coordinate->subjectNumber;

        // Get next coordinate;
        coordinate = index_getNextCoordinate();
    }

    // Perform second pass and extend hits found in index
    lastSequenceNumber = -1;
    sequenceCount = 0;
	descriptionStart = 0; descriptionLength = 0;
    databaseOffset = 0;

    coordinate = index_getFirstCoordinate();
    while (coordinate != NULL)
    {
		if (coordinate->subjectNumber != lastSequenceNumber)
        {
            address = startAddress - 40 + sequencePositions[coordinate->subjectNumber];
			descriptionStart = descriptionLocations[coordinate->subjectNumber];

            vbyte_getVbyte(address, &descriptionLength);
            vbyte_getVbyte(address, &subjectLength);

//            printf("%d) SeqPos=%d DescPos=%d Subject length=%d\n",
//                   coordinate->subjectNumber, index_sequencePositions[coordinate->subjectNumber],
//                   index_descriptionLocations[coordinate->subjectNumber], subjectLength);

            // Read a third vbyte; the total length of sequence data including wildcard info
            vbyte_getVbyte(address, &encodedLength);

            // Determine number of packed bytes and number of remaining letters
            numPackedBytes = subjectLength / 4;
            numRemaining = subjectLength % 4;

            subject = address;
            alignments_currentAlignment = NULL;
		}

//    	printf("Subject number=%d dloc=%d offset=%d,%d\n", coordinate->subjectNumber,
//               descriptionStart, coordinate->queryOffset, coordinate->subjectOffset);

        lastSequenceNumber = coordinate->subjectNumber;

        // Process the coordinate
        subjectPosition = coordinate->subjectOffset;
        queryPosition = coordinate->queryOffset;
        queryPosition -= subjectPosition % 4;
        subjectPosition /= 4;

        address = subject + subjectPosition;

        #ifdef VERBOSE
        if (parameters_verboseDloc == descriptionStart)
        	printf("Hit at %d,%d\n", queryPosition, subjectPosition);
		#endif

        #ifndef NO_STAGE2
        // Calculate the diagonal this hit is on
        diagonal = (subjectPosition * 4 - queryPosition + hitMatrix_queryLength)
                 & hitMatrix_diagonalMask;

        // If we have not extended past this point on this diagonal
        lastHit = hitMatrix_furthest + diagonal;

        if (*lastHit < address && ungappedExtension_checkHit(queryPosition, subjectPosition,
                PSSMatrix, subject, subjectLength, index_wordSize + index_intervalSize - 1))
        {
            #ifdef VERBOSE
            if (parameters_verboseDloc == descriptionStart)
        		printf("CheckHit score=%d\n", queryPosition, subjectPosition, ungappedExtension_bestScore);
            #endif

            // Perform ungapped extension
            queryPosition += 4; subjectPosition++;
            ungappedExtension
                = ungappedExtension_nucleotideExtend(queryPosition,
                  subjectPosition, PSSMatrix, subject, subjectLength);

            // Update furthest reached value for the diagonal
            *lastHit = ungappedExtension_subjectEndReached;

            // If extension scores high enough to trigger gapping
            if (ungappedExtension)
            {
//            	printf("Score=%d\n", ungappedExtension->nominalScore);
                // Increment count of number of trigger extensions
                blast_numTriggerExtensions++;

                // Create new alignment if needed
                if (alignments_currentAlignment == NULL)
                {
                    alignments_createNew(descriptionStart, descriptionLength, subject,
                                         subjectLength, encodedLength);
                }

                // Add this extension to the alignment
                alignments_addUngappedExtension(ungappedExtension);
            }

            blast_numUngappedExtensions++;
        }
        #endif

        // Get next coordinate;
        coordinate = index_getNextCoordinate();
    }

	free(sequencePositions);
	free(descriptionLocations);
}*/

// SSearch a protein database using Smith-Waterman algorithm
void search_proteinSsearch(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                           uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
	unsigned char *subject, *address;
	int4 subjectLength;
	struct gappedExtension* gappedExtension;
	struct dpResults dpResults, reverseDpResults;

    while (sequenceCount < numSequences)
    {
    	descriptionLength = sequenceData[sequenceCount].descriptionLength;
    	descriptionStart = sequenceData[sequenceCount].descriptionStart;
    	subjectLength = sequenceData[sequenceCount].sequenceLength;
    	encodedLength = sequenceData[sequenceCount].encodedLength;
		address = subject = sequenceData[sequenceCount].sequence;

        // New sequence, new possible alignment
        alignments_currentAlignment = NULL;

        // Perform SW score only
		dpResults = smithWatermanScoring_score(PSSMatrix, subjectLength, subject);

//        printf("%d\n", dpResults.bestScore);

        // If above e-value cutoff
		if (dpResults.bestScore >= blast_gappedNominalCutoff
            && alignments_isFinalAlignment(dpResults.bestScore))
		{
            // Perform SW alignment in the reverse direction, to find the start of the optimal alignment
        	reverseDpResults = smithWatermanScoring_scoreReverse(PSSMatrix, subjectLength, subject,
                                                                 dpResults.best);

            // Collect traceback information and store alignment
			gappedExtension = smithWatermanTraceback_build(PSSMatrix, subjectLength, subject,
                                                           reverseDpResults.best, dpResults.best);

            if (reverseDpResults.bestScore != dpResults.bestScore ||
                dpResults.bestScore != gappedExtension->nominalScore)
            {
                fprintf(stderr, "Error: Forward and reverse Smith-Waterman alignment scores do not match\n");
                exit(-1);
            }

            gappedExtension_score(gappedExtension);

			alignments_createNew(descriptionStart, descriptionLength, subject, subjectLength, encodedLength);
			alignments_addGappedExtension(alignments_currentAlignment, gappedExtension);
			alignments_addFinalAlignment(dpResults.bestScore, alignments_currentAlignment);
		}

        // Advance to next sequence
        address += encodedLength - 1;

        sequenceCount++;

        // Every so often print status .
        if ((sequenceCount % tickFrequency == 0) &&
            parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
        {
            #ifndef VERBOSE
            printf(".");
            fflush(stdout);
			#endif
        }
    }

	// Sort the final alignments by refined score
	alignments_sortFinalAlignments();
}

// Ssearch a nucleotide database using Smith-waterman algorithm
void search_nucleotideSsearch(struct PSSMatrix PSSMatrix, struct sequenceData* sequenceData,
                              uint4 numSequences, uint4 tickFrequency)
{
	uint4 sequenceCount = 0;
	uint4 descriptionStart = 0, descriptionLength = 0, encodedLength;
	unsigned char *subject, *unpackedSubject, *address;
	int4 subjectLength;
	struct gappedExtension* gappedExtension;
	struct dpResults dpResults, reverseDpResults;

    while (sequenceCount < numSequences)
    {
    	descriptionLength = sequenceData[sequenceCount].descriptionLength;
    	descriptionStart = sequenceData[sequenceCount].descriptionStart;
    	subjectLength = sequenceData[sequenceCount].sequenceLength;
    	encodedLength = sequenceData[sequenceCount].encodedLength;
		address = subject = sequenceData[sequenceCount].sequence;

        // New sequence, new possible alignment
        alignments_currentAlignment = NULL;

        // Unpack the subject
        unpackedSubject = encoding_byteUnpack(subject, subjectLength);

        // Re-insert wildcards
        encoding_insertWilds(unpackedSubject, subject + ((subjectLength + 3) / 4),
                             subject + encodedLength);

		// Perform SW score only
		dpResults = smithWatermanScoring_score(PSSMatrix, subjectLength, unpackedSubject);

        // If above e-value cutoff
		if (dpResults.bestScore >= blast_gappedNominalCutoff
            && alignments_isFinalAlignment(dpResults.bestScore))
		{
        	// Perform SW alignment in the reverse direction, to find the start of the optimal alignment
        	reverseDpResults = smithWatermanScoring_scoreReverse(PSSMatrix, subjectLength,
                               unpackedSubject, dpResults.best);

            // Collect traceback information and store alignment
			gappedExtension = smithWatermanTraceback_build(PSSMatrix, subjectLength, unpackedSubject,
                                                           reverseDpResults.best, dpResults.best);

            if (reverseDpResults.bestScore != dpResults.bestScore ||
                dpResults.bestScore != gappedExtension->nominalScore)
            {
                fprintf(stderr, "Error: Forward and reverse Smith-Waterman alignment scores do not match\n");
                fprintf(stderr, "Forward Score=%d End=%d,%d\n", dpResults.bestScore,
                        dpResults.best.queryOffset, dpResults.best.subjectOffset);
                fprintf(stderr, "Reverse Score=%d End=%d,%d\n", reverseDpResults.bestScore,
                        reverseDpResults.best.queryOffset, reverseDpResults.best.subjectOffset);
                fprintf(stderr, "Traceback Score=%d End=%d,%d\n", gappedExtension->nominalScore,
                        gappedExtension->queryEnd, gappedExtension->subjectEnd);
//                exit(-1);
            }

            gappedExtension_score(gappedExtension);

			alignments_createNew(descriptionStart, descriptionLength, unpackedSubject, subjectLength, 0);
			alignments_addGappedExtension(alignments_currentAlignment, gappedExtension);
			alignments_addFinalAlignment(dpResults.bestScore, alignments_currentAlignment);
		}
		else
        {
        	free(unpackedSubject);
        }

        sequenceCount++;

        // Every so often print status .
        if ((sequenceCount % tickFrequency == 0) &&
            parameters_outputType != parameters_xml && parameters_outputType != parameters_tabular)
        {
            #ifndef VERBOSE
            printf(".");
            fflush(stdout);
            #endif
        }
    }
}

/*
// Testing code
void search_compareScoring(struct alignment* alignment, struct PSSMatrix PSSMatrix)
{
	struct ungappedExtension* ungappedExtension;
	int4 score1, score2, score3, score4, best1 = 0, best3 = 0, count;

	// For each ungapped extension for this subject, highest-scoring first
	ungappedExtension = alignment->ungappedExtensions;
	while (ungappedExtension != NULL)
	{
		blast_dloc = alignment->descriptionLocation;

    	ungappedExtension_findSeed(ungappedExtension, PSSMatrix, alignment->subject);

        count = 0;
        while (count < 5)
        {
        score1
        	= tableGappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
		                               alignment->subject, statistics_gappedNominalDropoff);
		count++;
        }

        score2
			= nuGappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
		                            alignment->subject, statistics_gappedNominalDropoff);

        score4
        	= bytepackGappedScoring_score(ungappedExtension, PSSMatrix, alignment->subjectLength,
		        alignment->subject, statistics_gappedNominalDropoff);


		printf("UE=%d TableDriven=%d GappedScoring=%d Bytepacked=%d Dloc=%d Slength=%d SeedQ=%d SeedS=%d\n",
			ungappedExtension->nominalScore, score1, score2, score4,
            alignment->descriptionLocation, alignment->subjectLength,
			ungappedExtension->seed.queryOffset, ungappedExtension->seed.subjectOffset);

        // Remove poor-scoring ungapped extension from start of list, and move
		// to the now-first ungapped extension for processing
		alignment->ungappedExtensions = ungappedExtension->next;
		ungappedExtension = alignment->ungappedExtensions;
	}
}

void search_compareScorings(struct PSSMatrix PSSMatrix)
{
	struct alignment* alignment;

    memBlocks_resetCurrent(alignments_alignments);

	while ((alignment = (struct alignment*)memBlocks_getCurrent(alignments_alignments)))
	{
        search_compareScoring(alignment, PSSMatrix);
	}
}
*/

