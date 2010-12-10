extern unsigned char * wordLookupDFA;
extern struct groupFP *wordLookupDFA_groupsFP;
TIMERECORD timeRecord;

//texture<unsigned char, 1> texSubjectSequences;

__global struct parameters {
	char 	parameters_wordSize;
	unsigned char encoding_numCodes;
	char    parameters_overlap;
	int4 	wordLookupDFA_numCodes;
	uint4	additionalQueryPositionOffset;
	int4	statistics_ungappedNominalDropoff;
	int4	blast_ungappedNominalTrigger;
	int4    parameters_A;
	uint4	ungappedExtensionsPerThread;
	uint4	ungappedExtAdditionalStartLoc;
};

#define TARGET_THREAD 0
#define UNGAPEXT_PER_THREAD 150
#define TOTAL_UNGAPPED_EXT 1500000

//NOOOOOOOOOO idea what to do about these
//Looks like this first part is all Kernel Stuff, will need forking out
__global struct ungappedExtension *ungappedExtension_oneHitExtendD(unsigned char *, int4, unsigned
char*, 
	struct PSSMatrixFP, unsigned char*, uint4 *, unsigned char, int4, int4,
	struct ungappedExtension *, uint4 *, uint4, int4);
struct coordinate ungappedExtension_findProteinSeed(struct ungappedExtension*,
		   struct PSSMatrixFP, unsigned char*, unsigned char);

//__device__ __constant__ struct groupFP wordLookupDFA_groupsC[400];
//This doesn't seem to be used right now oddly enough
//__constant const int2 scoreMatrixC[1640];
__constant unsigned char querySequenceC[40000];
uint4 global_sequenceCount;
uint4 global_numAdditionalTriggerExtensions;

//Embarrassingly parallel approach is used. One thread is used for
//the hit detection of one sequence
__kernel void search_protein1hitKernel(struct PSSMatrixFP *PSSMatrixFP,
										 int2 *matrixBody,
									  	 struct sequenceDataFP *sequenceDataFP, 
										 unsigned char *sequence,
										 struct parameters *parametersFP,
										 struct groupFP *wordLookupDFA_groupsFP,
										 unsigned char *wordLookupDFAFP,
										 uint4 *blast_numUngappedExtensions,
										 uint4 *blast_numTriggerExtensions,
										 uint4 *blast_numHits,
										 uint4 *hitMatrix_furthestp,
										 uint4 *hitMatrix_offsetp,
										 struct ungappedExtension * ungappedExtension_extensionsp,
										 uint4 nTotalSequenceNum)
{
	int bid = get_local_id(1) * get_global_size(2) + get_local_id(2);
//	int bid = blockIdx.x * gridDim.y + blockIdx.y;
	int tid = bid * get_local_size(1) * get_local_size(2) + get_group_id(1) * get_local_size(2) + get_group_id(2);
//	int tid = bid * blockDim.x * blockDim.y + threadIdx.x * blockDim.y + threadIdx.y;

	unsigned char *subject, *sequenceEnd, *address;
	int4 subjectOffset, count;
	unsigned char currentWord, *currentBlock;
	struct groupFP *currentGroupFP;
	uint2 *wordLookupDFA_AddiPositions;
	uint4 numOfTriggerExtensions = 0;
	uint2 *queryOffsets, queryOffset;
	struct ungappedExtension* ungappedExtension_current;
	int4 diagonal;
	uint4 *lastHitFP;
	uint4 ungappedExtension_subjectEndReachedFP;
	uint4 *hitMatrix_Local;
	uint4 sequenceCount;

	hitMatrix_Local = hitMatrix_furthestp +  hitMatrix_offsetp[tid] + PSSMatrixFP->length;
	ungappedExtension_extensionsp->start.subjectOffset = 0;
	ungappedExtension_current = ungappedExtension_extensionsp + tid * UNGAPEXT_PER_THREAD;
	wordLookupDFA_AddiPositions = (uint2 *)((char *)wordLookupDFAFP + 
								  parametersFP->additionalQueryPositionOffset);

	//Set the PSSMatrix body
	PSSMatrixFP->matrix = matrixBody + parametersFP->encoding_numCodes;

	sequenceCount = tid;
	while (sequenceCount < nTotalSequenceNum)
	{
		subject = address = sequence + sequenceDataFP[sequenceCount].offset;
		
		if (sequenceDataFP[sequenceCount].sequenceLength >= parametersFP->parameters_wordSize)
		{
			currentGroupFP = wordLookupDFA_groupsFP;
			//currentGroupFP = wordLookupDFA_groupsC;

			count = 1;
			while (count < parametersFP->parameters_wordSize)
			{
				if (*address < parametersFP->wordLookupDFA_numCodes)
				{
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
					//currentGroupFP = &wordLookupDFA_groupsC[currentGroupFP->nextGroups + *address];
				}
				else
				{
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
					//currentGroupFP = &wordLookupDFA_groupsC[currentGroupFP->nextGroups];
				}
				address++;
				count++;
			}

			sequenceEnd = subject + sequenceDataFP[sequenceCount].sequenceLength;

			while (address < sequenceEnd)
			{
				currentBlock = &wordLookupDFAFP[currentGroupFP->nextWords];

				// If current code is a regular letter
				if (*address < parametersFP->wordLookupDFA_numCodes)
				{
					
					currentWord = currentBlock[*address];
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
					//currentGroupFP = &wordLookupDFA_groupsC[currentGroupFP->nextGroups + *address];
				}
				else
				{
					if (address >= sequenceEnd)
						break;

					currentWord = *currentBlock;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
					//currentGroupFP = &wordLookupDFA_groupsC[currentGroupFP->nextGroups];
				}
				
				if (currentWord)
				{
					subjectOffset = address - subject;
					// At least one query position, stored at an extenal address
					queryOffsets = ((uint2*)currentBlock) - currentWord;
					
					if (!*queryOffsets)
					{
						// Go to an outside address for additional positions
						queryOffsets = wordLookupDFA_AddiPositions
									+ (*(queryOffsets + 1) * constants_max_int2) + *(queryOffsets + 2);
					}

					do
					{
						queryOffset = *queryOffsets;
		
						#ifndef NO_STAGE2
						// Calculate the diagonal this hit is on
						diagonal = subjectOffset - queryOffset;

						// If we have not extended past this point on this diagonal
						lastHitFP = hitMatrix_Local + diagonal;
						
						if (*lastHitFP < address - sequence)
						{
							//Number of extensions for each subject sequence
							blast_numUngappedExtensions[tid] ++;

							// If only one hit triggered this extension
							ungappedExtension_oneHitExtendD(sequence, 
									queryOffset,
									address, *PSSMatrixFP, subject, 
									&ungappedExtension_subjectEndReachedFP,
									parametersFP->encoding_numCodes,
									parametersFP->statistics_ungappedNominalDropoff,
									parametersFP->blast_ungappedNominalTrigger,
									ungappedExtension_current,
									&numOfTriggerExtensions,
									sequenceCount,
									tid);

							// Update furthest reached value for the diagonal
							*lastHitFP = ungappedExtension_subjectEndReachedFP;
						}
						#endif
					
						queryOffsets++; blast_numHits[tid]++;
					}
					while(*queryOffsets);
				}
				address++;
			}
		}
		//option=======================================================
		//sequenceCount = atomicAdd(&global_sequenceCount, 1);
		sequenceCount += get_global_size(1) * get_local_size(1);
//		sequenceCount += gridDim.x * blockDim.x;
		//============================================================
	}

	blast_numTriggerExtensions[tid] = (uint4) numOfTriggerExtensions;
	return;
}

__global struct ungappedExtension *ungappedExtension_oneHitExtendD(unsigned char* sequenceStart,
						int4 queryOffset, 
						unsigned char *subjectHit, 
						struct PSSMatrixFP PSSMatrixFP, unsigned char *subject,
						uint4 *sequenceHitEnd, unsigned char encoding_numCodes,
						int4 statistics_ungappedNominalDropoff,
						int4 blast_ungappedNominalTrigger,
						struct ungappedExtension *ungappedExtension_extensions,
						uint4 *numOfTriggerExtensions,
						uint4 sequenceCount,
						int4 tid)
{
	int2* queryPosition;
	//int4 queryPosition;
	unsigned char* subjectPosition, *subjectStart, *subjectEnd;
	int4 changeSinceBest = 0;
	int4 dropoff, originalDropoff;
	int4 ungappedExtension_bestScore;

	originalDropoff = dropoff = -statistics_ungappedNominalDropoff;	
	ungappedExtension_bestScore = 0;

	// Start at queryEnd,subjectEnd (right/last hit position)
	queryPosition = PSSMatrixFP.matrix + queryOffset * encoding_numCodes;
	//queryPosition = queryOffset + 1;
	subjectPosition = subjectStart = subjectHit;
	
	while (changeSinceBest > dropoff)
	{
		changeSinceBest += queryPosition[*subjectPosition];
		//changeSinceBest += scoreMatrixC[querySequenceC[queryPosition] * encoding_numCodes + (*subjectPosition)];

		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				queryPosition = queryPosition - encoding_numCodes;
				//queryPosition = queryPosition - 1;
				subjectPosition--;

				changeSinceBest = queryPosition[*subjectPosition];
				//changeSinceBest = scoreMatrixC[querySequenceC[queryPosition] * encoding_numCodes + (*subjectPosition)];
			}
			while (changeSinceBest > 0);

			subjectStart = subjectPosition;
		}

		queryPosition = queryPosition - encoding_numCodes;
		//queryPosition = queryPosition - 1;
		subjectPosition--;
	}
	
	// Correct for extra decrement
	subjectStart++;

	// Starting at right/last hit position again
	queryPosition = PSSMatrixFP.matrix + (queryOffset + 1) * encoding_numCodes;
	//queryPosition = (queryOffset + 2);
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
		//changeSinceBest += scoreMatrixC[querySequenceC[queryPosition] * encoding_numCodes + (*subjectPosition)];
		
		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				queryPosition = queryPosition + encoding_numCodes;
				//queryPosition = queryPosition + 1;

				subjectPosition++;
				changeSinceBest = queryPosition[*subjectPosition];
				//changeSinceBest = scoreMatrixC[querySequenceC[queryPosition] * encoding_numCodes + (*subjectPosition)];
			}
			while (changeSinceBest > 0);

			subjectEnd = subjectPosition;
			
			// Check need for change in dropoff
			if ((dropoff = -ungappedExtension_bestScore) < originalDropoff)
			{
				dropoff = originalDropoff;
			}
		}
		queryPosition = queryPosition + encoding_numCodes;
		//queryPosition = queryPosition + 1;

		subjectPosition++;
	}
	
	subjectEnd--;

	//*sequenceHitEnd = subjectPosition - subject;
	*sequenceHitEnd = subjectPosition - sequenceStart;

	if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger)
	{
		int4 diagonal;
		struct ungappedExtension* newUngappedExtension = NULL;
		
		newUngappedExtension = &ungappedExtension_extensions[*numOfTriggerExtensions];
		// Calculate diagonal
		diagonal = (subjectHit - subject) - queryOffset;

		// Determine offsets from pointers
		newUngappedExtension->start.subjectOffset = subjectStart - subject;
		newUngappedExtension->end.subjectOffset = subjectEnd - subject;
		newUngappedExtension->start.queryOffset = newUngappedExtension->start.subjectOffset - diagonal;
		newUngappedExtension->end.queryOffset = newUngappedExtension->end.subjectOffset - diagonal;
		newUngappedExtension->seed = ungappedExtension_findProteinSeed(newUngappedExtension, 
									 PSSMatrixFP, subject, encoding_numCodes);

		newUngappedExtension->next = NULL;
		newUngappedExtension->nominalScore = ungappedExtension_bestScore;
		newUngappedExtension->status = ungappedExtension_UNGAPPED;
		newUngappedExtension->sequenceCount = sequenceCount;

		//Shucai
		//Record the number of hits satisfying the next step
		(*numOfTriggerExtensions)++;

		return newUngappedExtension;
	}
	else
	{
		return NULL;
	}

}

__global struct coordinate ungappedExtension_findProteinSeed(
					struct ungappedExtension* ungappedExtension,
					struct PSSMatrixFP PSSMatrixFP, 
					unsigned char* subject, 
					unsigned char encoding_numCodes)
{
	int2 *queryWindowStart, *queryWindowEnd;
	unsigned char *subjectWindowStart, *subjectWindowEnd;

	int2* bestQueryPosition;
	unsigned char* bestSubjectPosition;
	int4 bestSegmentScore;
	int4 nominalScore, count;
	struct coordinate seed;

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
		queryWindowStart = queryWindowEnd = PSSMatrixFP.matrix + ungappedExtension->start.queryOffset * encoding_numCodes;
		subjectWindowStart = subjectWindowEnd = subject + ungappedExtension->start.subjectOffset;

		// Find initial score for first 11 positions
		nominalScore = 0;
		count = 0;
		while (count < 11)
		{
			nominalScore += queryWindowEnd[*subjectWindowEnd];
			queryWindowEnd += encoding_numCodes;
			subjectWindowEnd++;
			count++;
		}

		queryWindowEnd -= encoding_numCodes;
		subjectWindowEnd--;

		// By default first-11 positions gives best position and score
		bestQueryPosition = queryWindowStart;
		bestSubjectPosition = subjectWindowStart;
		bestSegmentScore = nominalScore;

		// Now slide the window across and record the better scores/positions
		while (queryWindowEnd < PSSMatrixFP.matrix + ungappedExtension->end.queryOffset * encoding_numCodes)
		{
			// Advance window end, add new position value
			queryWindowEnd += encoding_numCodes;
			subjectWindowEnd++;

			nominalScore += queryWindowEnd[*subjectWindowEnd];
			// Remove position that we will leave behind
			nominalScore -= queryWindowStart[*subjectWindowStart];

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
		seed.queryOffset = (bestQueryPosition - PSSMatrixFP.matrix) / encoding_numCodes + 5;
		seed.subjectOffset = bestSubjectPosition + 5 - subject;
	}

	return seed;
}

/*************************************************************************
 *  The following functions are for 2-hit extension
 *************************************************************************/
__global struct ungappedExtension *ungappedExtension_twoHitExtendD(unsigned char* sequenceStart,
						int4 queryOffset, 
						unsigned char *subjectHit,
						uint4 lastHitFP,
						struct PSSMatrixFP PSSMatrixFP, unsigned char *subject,
						uint4 *sequenceHitEnd, unsigned char encoding_numCodes,
						int4 statistics_ungappedNominalDropoff,
						int4 blast_ungappedNominalTrigger,
						int4 ungappedExtensionsPerThread,
						struct ungappedExtension *ungappedExtension_extensions,
						struct ungappedExtension *ungappedExtension_additonal,
						uint4 *numOfTriggerExtensions,
						uint4 sequenceCount,
						int4 tid)
{
	int2* queryPosition;
	unsigned char* subjectPosition, *subjectStart, *subjectEnd;
	int4 changeSinceBest = 0;
	int4 dropoff, originalDropoff;
	int4 ungappedExtension_bestScore;

	originalDropoff = dropoff = -statistics_ungappedNominalDropoff;	
	ungappedExtension_bestScore = 0;

	// Start at queryEnd,subjectEnd (right/last hit position)
	queryPosition = PSSMatrixFP.matrix + queryOffset * encoding_numCodes;
	subjectPosition = subjectStart = subjectHit;
	
	while (changeSinceBest > dropoff)
	{
		changeSinceBest += queryPosition[*subjectPosition];

		// If we have got a positive score
		if (changeSinceBest > 0)
		{
			// Keep updating best score and resetting change-since-best
			// whilst we are reading positive scores
			do
			{
				ungappedExtension_bestScore += changeSinceBest;
				queryPosition = queryPosition - encoding_numCodes;
				subjectPosition--;

				changeSinceBest = queryPosition[*subjectPosition];
			}
			while (changeSinceBest > 0);

			subjectStart = subjectPosition;
		}

		queryPosition = queryPosition - encoding_numCodes;
		subjectPosition--;
	}

	// Correct for extra decrement
	subjectStart++;

	if (subjectStart - sequenceStart > lastHitFP)
	{
		*sequenceHitEnd = subjectHit - sequenceStart;
		return NULL;
	}

	// Starting at right/last hit position again
	queryPosition = PSSMatrixFP.matrix + (queryOffset + 1) * encoding_numCodes;
	subjectPosition = subjectHit + 1;
	subjectEnd = subjectHit;
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
				queryPosition = queryPosition + encoding_numCodes;
				subjectPosition++;
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
		queryPosition = queryPosition + encoding_numCodes;
		subjectPosition++;
	}

	subjectEnd--;

	*sequenceHitEnd = subjectEnd - sequenceStart;

	if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger)
	{
		int4 diagonal;
		struct ungappedExtension* newUngappedExtension = NULL;
		
		if (*numOfTriggerExtensions >= ungappedExtensionsPerThread)
		{
			newUngappedExtension = &ungappedExtension_additonal[atomicAdd(&global_numAdditionalTriggerExtensions, 1)];
		}
		else
		{
			newUngappedExtension = &ungappedExtension_extensions[*numOfTriggerExtensions];
		}

		// Calculate diagonal
		diagonal = (subjectHit - subject) - queryOffset;

		// Determine offsets from pointers
		newUngappedExtension->start.subjectOffset = subjectStart - subject;
		newUngappedExtension->end.subjectOffset = subjectEnd - subject;
		newUngappedExtension->start.queryOffset = newUngappedExtension->start.subjectOffset - diagonal;
		newUngappedExtension->end.queryOffset = newUngappedExtension->end.subjectOffset - diagonal;
		newUngappedExtension->seed = ungappedExtension_findProteinSeed(newUngappedExtension, 
									 PSSMatrixFP, subject, encoding_numCodes);

		newUngappedExtension->next = NULL;
		newUngappedExtension->nominalScore = ungappedExtension_bestScore;
		newUngappedExtension->status = ungappedExtension_UNGAPPED;
		newUngappedExtension->sequenceCount = sequenceCount;
		newUngappedExtension->tid = tid;

		//Shucai
		//Record the number of hits satisfying the next step
		(*numOfTriggerExtensions)++;

		return newUngappedExtension;
	}
	else
	{
		return NULL;
	}

}

//Embarrassingly parallel approach is used. One thread is used for
//the hit detection of one sequence
__kernel void search_protein2hitKernel(struct PSSMatrixFP *PSSMatrixFP,
										 int2 *matrixBody,
									  	 struct sequenceDataFP *sequenceDataFP, 
										 unsigned char *sequence,
										 struct parameters *parametersFP,
										 struct groupFP *wordLookupDFA_groupsFP,
										 unsigned char *wordLookupDFAFP,
										 uint4 *blast_numUngappedExtensions,
										 uint4 *blast_numTriggerExtensions,
										 uint4 *blast_numHits,
										 uint4 *hitMatrix_furthestp,
										 uint4 *hitMatrix_offsetp,
										 struct ungappedExtension * ungappedExtension_extensionsp,
										 uint4 nTotalSequenceNum)
{
	int tid = get_local_id(1) * get_local_size(1) + get_group_id(x);
//	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	unsigned char *subject, *sequenceEnd, *address;
	int4 subjectOffset, count;
	unsigned char currentWord, *currentBlock;
	struct groupFP *currentGroupFP;
	uint2 *wordLookupDFA_AddiPositions;
	uint4 numOfTriggerExtensions = 0;
	uint2 *queryOffsets, queryOffset;
	struct ungappedExtension* ungappedExtension_current;
	struct ungappedExtension* ungappedExtension_additional;
	int4 diagonal;
	uint4 *lastHitFP;
	uint4 ungappedExtension_subjectEndReachedFP;
	uint4 *hitMatrix_Local;
	uint4 sequenceCount;
	int4  distance;

	hitMatrix_Local = hitMatrix_furthestp +  hitMatrix_offsetp[tid] + PSSMatrixFP->length;
	//ungappedExtension_current = ungappedExtension_extensionsp + tid * UNGAPEXT_PER_THREAD;
	ungappedExtension_extensionsp->start.subjectOffset = 0;
	ungappedExtension_current = ungappedExtension_extensionsp + tid * parametersFP->ungappedExtensionsPerThread;
	ungappedExtension_additional = ungappedExtension_extensionsp + parametersFP->ungappedExtAdditionalStartLoc;
	wordLookupDFA_AddiPositions = (uint2 *)((char *)wordLookupDFAFP + 
								  parametersFP->additionalQueryPositionOffset);

	//Set the PSSMatrix body
	PSSMatrixFP->matrix = matrixBody + parametersFP->encoding_numCodes;

	sequenceCount = tid;
	while (sequenceCount < nTotalSequenceNum)
	{
		subject = address = sequence + sequenceDataFP[sequenceCount].offset;
		
		if (sequenceDataFP[sequenceCount].sequenceLength >= parametersFP->parameters_wordSize)
		{
			currentGroupFP = wordLookupDFA_groupsFP;

			count = 1;
			while (count < parametersFP->parameters_wordSize)
			{
				if (*address < parametersFP->wordLookupDFA_numCodes)
				{
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
				}
				else
				{
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
				}
				address++;
				count++;
			}

			sequenceEnd = subject + sequenceDataFP[sequenceCount].sequenceLength;

			while (address < sequenceEnd)
			{
				currentBlock = &wordLookupDFAFP[currentGroupFP->nextWords];

				// If current code is a regular letter
				if (*address < parametersFP->wordLookupDFA_numCodes)
				{
					currentWord = currentBlock[*address];
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
				}
				else
				{
					if (address >= sequenceEnd)
						break;

					currentWord = *currentBlock;
					currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
				}
				
				if (currentWord)
				{
					subjectOffset = address - subject;
					// At least one query position, stored at an extenal address
					queryOffsets = ((uint2*)currentBlock) - currentWord;
					
					if (!*queryOffsets)
					{
						// Go to an outside address for additional positions
						queryOffsets = wordLookupDFA_AddiPositions
									+ (*(queryOffsets + 1) * constants_max_int2) + *(queryOffsets + 2);
					}

					do
					{
						queryOffset = *queryOffsets;
		
						#ifndef NO_STAGE2
						// Calculate the diagonal this hit is on
						diagonal = subjectOffset - queryOffset;

						// If we have not extended past this point on this diagonal
						lastHitFP = hitMatrix_Local + diagonal;

						distance = (address - sequence) - *lastHitFP;
						if (distance >= parametersFP->parameters_A)
						{
							*lastHitFP = address - sequence;
						}
						else if (distance >= parametersFP->parameters_overlap)
						{
							//Number of extensions for each subject sequence
							blast_numUngappedExtensions[tid] ++;

							// If only one hit triggered this extension
							ungappedExtension_twoHitExtendD(sequence, 
									queryOffset,
									address, *lastHitFP, *PSSMatrixFP, subject, 
									&ungappedExtension_subjectEndReachedFP,
									parametersFP->encoding_numCodes,
									parametersFP->statistics_ungappedNominalDropoff,
									parametersFP->blast_ungappedNominalTrigger,
									parametersFP->ungappedExtensionsPerThread,
									ungappedExtension_current,
									ungappedExtension_additional,
									&numOfTriggerExtensions,
									sequenceCount,
									tid);

							// Update furthest reached value for the diagonal
							*lastHitFP = ungappedExtension_subjectEndReachedFP;
						}
						#endif
					
						queryOffsets++; blast_numHits[tid]++;
					}
					while(*queryOffsets);
				}
				address++;
			}
		}

		//option=======================================================
		//sequenceCount = atomicAdd(&global_sequenceCount, 1);
		sequenceCount += get_global_size(1) * get_local_size(1);
//		sequenceCount += gridDim.x * blockDim.x;
		//============================================================
	}

	blast_numTriggerExtensions[tid] = (uint4) numOfTriggerExtensions;
	return;
}
