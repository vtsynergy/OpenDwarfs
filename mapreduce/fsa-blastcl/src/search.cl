//#include "../include/global.h"
//#include "../include/constants.h"

#pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_amd_printf : enable
#define ungappedExtension_UNGAPPED 1
#define NULL 0
#define constants_max_short 32767
#define TARGET_THREAD 0
#define UNGAPEXT_PER_THREAD 150
#define TOTAL_UNGAPPED_EXT 1500000
struct __attribute__((aligned (8))) PSSMatrixFP
{
	int length;
	int strandLength;
	int highestValue;
	int lowestValue;
__global	short* matrix;
    unsigned char* queryCodes;
	unsigned char* bestMatchCodes;
	unsigned char* bytePackedCodes;
	unsigned char* xorCodes;
};


struct __attribute__((aligned (8))) sequenceDataFP
{
	uint sequenceLength;
	uint descriptionStart;
	uint descriptionLength;
	uint encodedLength;
    uint offset;
};

typedef struct __attribute__((aligned (8))) strTimeRecord {
	uint iniTime;
	uint preProcessTime;
	uint dataCopyTimeH2D;
	uint searchTime;
	uint dataCopyTimeD2H;
	uint addUngappedExtensionTime;
	uint postProcessTime;
	uint hitUngappedExtTime;
	uint gappedAlignmentTime;
	uint finalAlignmentTime;
	uint totalTime;
} TIMERECORD;

struct __attribute__((aligned (8))) coordinate
{
    int queryOffset;
	int subjectOffset;
};

struct __attribute__ ((aligned (8))) ungappedExtension
{
	int nominalScore;
	//Shucai
	uint sequenceCount;
	uint tid;
	struct ungappedExtension* next;
	struct coordinate start;
	struct coordinate end;
	struct coordinate seed;
	int status; //char status;
//	char pad[11];
};

//__global unsigned char* wordLookupDFA;
//__constant unsigned char * wordLookupDFA;
struct __attribute__ ((aligned (8))) groupFP //for parallelization
{
	int nextWords;
	int nextGroups;
};
//__constant struct groupFP  *wordLookupDFA_groupsFP;
__constant TIMERECORD timeRecord;

//texture<unsigned char, 1> texSubjectSequences;

struct __attribute__ ((aligned (8))) parameters {
	int 	wordLookupDFA_numCodes;
	uint	additionalQueryPositionOffset;
	int	statistics_ungappedNominalDropoff;
	int	blast_ungappedNominalTrigger;
	int    parameters_A;
	uint	ungappedExtensionsPerThread;
	uint	ungappedExtAdditionalStartLoc;
	char 	parameters_wordSize;
	unsigned char encoding_numCodes;
	char    parameters_overlap;
};


__global struct ungappedExtension *ungappedExtension_oneHitExtendD(__global unsigned char *, int, __global unsigned
char*, 
	struct PSSMatrixFP, __global unsigned char*, uint *, unsigned char, int, int,
	__global struct ungappedExtension *, uint *, uint, int);
struct coordinate ungappedExtension_findProteinSeed(__global struct ungappedExtension*,
		   struct PSSMatrixFP, __global unsigned char*, unsigned char);

//__device__ __constant__ struct groupFP wordLookupDFA_groupsC[400];
//This doesn't seem to be used right now oddly enough
//__constant const short scoreMatrixC[1640];
//__constant unsigned char querySequenceC[40000];
//__constant uint global_sequenceCount;
//__global uint global_numAdditionalTriggerExtensions;

__kernel void memSet(uint value, __global uint *mem) {
    mem[get_global_id(0)] = value;
//if (get_global_id(0) == 0) printf("MEMSET CHECK TID0 ADDRESS: %x %d\n", &mem[0], value);
}


/*************************************************************************
 *  The following functions are for 2-hit extension
 *************************************************************************/
__global struct ungappedExtension* ungappedExtension_twoHitExtendD(__global unsigned char* sequenceStart,
						int queryOffset, 
						__global unsigned char *subjectHit,
						uint lastHitFP,
						struct PSSMatrixFP PSSMatrixFP, __global unsigned char *subject,
						uint *sequenceHitEnd, unsigned char encoding_numCodes,
						int statistics_ungappedNominalDropoff,
						int blast_ungappedNominalTrigger,
						int ungappedExtensionsPerThread,
						__global struct ungappedExtension *ungappedExtension_extensions,
						__global struct ungappedExtension *ungappedExtension_additonal,
						uint *numOfTriggerExtensions,
						uint sequenceCount,
						int tid,
						__global uint* global_numAdditionalTriggerExtensions)
{
//printf("mrower\n");
//printf("UNGAPPED KERNEL %d\n", sizeof(struct ungappedExtension));

__global	short* queryPosition;
	__global unsigned char* subjectPosition, *subjectStart, *subjectEnd;
	int changeSinceBest = 0;
	int dropoff, originalDropoff;
	int ungappedExtension_bestScore;

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
//		printf("bower tid:%d %d %d\n", tid, ungappedExtension_bestScore, blast_ungappedNominalTrigger);
		int diagonal;
//		struct ungappedExtension* newUngappedExtension = 0;
__global struct ungappedExtension* newUngappedExtension = NULL;
		
		if ((*numOfTriggerExtensions) >= ungappedExtensionsPerThread)
		{
//			newUngappedExtension = &ungappedExtension_additonal[atomic_add(global_numAdditionalTriggerExtensions, 1)];
			newUngappedExtension = &ungappedExtension_additonal[atom_add(global_numAdditionalTriggerExtensions, 1)]; //Nvidia needs this line since they still only support OpenCL 1.0
		}
		else
		{
			newUngappedExtension = &ungappedExtension_extensions[(*numOfTriggerExtensions)];
		}

		// Calculate diagonal
		diagonal = (subjectHit - subject) - queryOffset;

//printf("diagonalu %d\n", diagonal);
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
//printf("scores %d %d %d %d\n", sequenceCount, ungappedExtension_bestScore, blast_ungappedNominalTrigger, newUngappedExtension->status);
		newUngappedExtension->tid = tid;

		//Shucai
		//Record the number of hits satisfying the next step
		(*numOfTriggerExtensions)++;
		

		return newUngappedExtension;
	}
	else
	{
		return 0;
//		return NULL;
	}

}

//Embarrassingly parallel approach is used. One thread is used for
//the hit detection of one sequence
__kernel void search_protein1hitKernel(__global struct PSSMatrixFP *PSSMatrixFP,
										 __global short *matrixBody,
									  	 __global struct sequenceDataFP *sequenceDataFP, 
										 __global unsigned char *sequence,
										 __global struct parameters *parametersFP,
										 __global struct groupFP *wordLookupDFA_groupsFP,
										 __global unsigned char *wordLookupDFAFP,
										 __global uint *blast_numUngappedExtensions,
										 __global uint *blast_numTriggerExtensions,
										 __global uint *blast_numHits,
										 __global uint *hitMatrix_furthestp,
										 __global uint *hitMatrix_offsetp,
										 __global struct ungappedExtension * ungappedExtension_extensionsp,
										 uint nTotalSequenceNum)
{
	int bid = get_local_id(0) * get_global_size(1) + get_local_id(1);
//	int bid = blockIdx * gridDim.y + blockIdx.y;
	int tid = bid * get_local_size(0) * get_local_size(1) + get_group_id(0) * get_local_size(1) + get_group_id(1);
//	int tid = bid * blockDim * blockDim.y + threadIdx * blockDim.y + threadIdx.y;

	__global unsigned char *subject, *sequenceEnd, *address;
	int subjectOffset, count;
unsigned char currentWord;
	__global unsigned char *currentBlock;
	__global struct groupFP *currentGroupFP;
	__global ushort *wordLookupDFA_AddiPositions;
	uint numOfTriggerExtensions = 0;
	__global ushort *queryOffsets;
	ushort queryOffset;
	__global struct ungappedExtension* ungappedExtension_current;
	int diagonal;
	__global uint *lastHitFP;
	uint ungappedExtension_subjectEndReachedFP;
	__global uint *hitMatrix_Local;
	uint sequenceCount;


	hitMatrix_Local = hitMatrix_furthestp + hitMatrix_offsetp[tid] + PSSMatrixFP->length;
	ungappedExtension_extensionsp->start.subjectOffset = 0;
	ungappedExtension_current = ungappedExtension_extensionsp + tid * UNGAPEXT_PER_THREAD;
	wordLookupDFA_AddiPositions = (__global unsigned short *)((__global char *)wordLookupDFAFP + parametersFP->additionalQueryPositionOffset);

	//Set the PSSMatrix body
PSSMatrixFP->matrix= matrixBody + parametersFP->encoding_numCodes;//	PSSMatrixFP->matrix = (short *)((int)(matrixBody + parametersFP->encoding_numCodes));

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
					queryOffsets = ((__global ushort*)currentBlock) - currentWord;
					
					if (!(*queryOffsets))
					{
						// Go to an outside address for additional positions
*queryOffsets = *wordLookupDFA_AddiPositions
                                     + (*(queryOffsets + 1) * constants_max_short) + *(queryOffsets + 2);
//*queryOffsets = *wordLookupDFA_AddiPositions
//                                     + (queryOffsets[1] * constants_max_short) + queryOffsets[2];
						//queryOffsets = wordLookupDFA_AddiPositions + (*((ushort*)((int)queryOffsets + 1 * sizeof(ushort*))) * constants_max_short) + (*(queryOffsets + 2));
					}

					do
					{
						queryOffset = *queryOffsets;
		
						#ifndef NO_STAGE2
						// Calculate the diagonal this hit is on
						diagonal = subjectOffset - queryOffset;

						// If we have not extended past this point on this diagonal
						lastHitFP = hitMatrix_Local + diagonal;
						
						if ((*lastHitFP) < address - sequence)
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
							*lastHitFP = ungappedExtension_subjectEndReachedFP; //VUVE
						}
						#endif
					
						queryOffsets++; blast_numHits[tid]++;
					}
					while((*queryOffsets));
				}
				address++;
			}
		}
		//option=======================================================
		//sequenceCount = atomicAdd(&global_sequenceCount, 1);
		sequenceCount += get_num_groups(0) * get_local_size(0);
//		sequenceCount += gridDim * blockDim;
		//============================================================
	}

	blast_numTriggerExtensions[tid] = (uint) numOfTriggerExtensions;
	return;
}

__global struct ungappedExtension *ungappedExtension_oneHitExtendD(__global unsigned char* sequenceStart,
						int queryOffset, 
						__global unsigned char *subjectHit, 
						struct PSSMatrixFP PSSMatrixFP, __global unsigned char *subject,
						uint *sequenceHitEnd, unsigned char encoding_numCodes,
						int statistics_ungappedNominalDropoff,
						int blast_ungappedNominalTrigger,
						__global struct ungappedExtension *ungappedExtension_extensions,
						uint *numOfTriggerExtensions,
						uint sequenceCount,
						int tid)
{
	__global short* queryPosition;
	//int queryPosition;
	__global unsigned char* subjectPosition, *subjectStart, *subjectEnd;
	int changeSinceBest = 0;
	int dropoff, originalDropoff;
	int ungappedExtension_bestScore;

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
		int diagonal;
		__global struct ungappedExtension* newUngappedExtension = NULL;
		
		newUngappedExtension = &ungappedExtension_extensions[(*numOfTriggerExtensions)];
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

struct coordinate ungappedExtension_findProteinSeed(
					__global struct ungappedExtension* ungappedExtension,
					struct PSSMatrixFP PSSMatrixFP, 
					__global unsigned char* subject, 
					unsigned char encoding_numCodes)
{
	__global short *queryWindowStart, *queryWindowEnd;
	__global unsigned char *subjectWindowStart, *subjectWindowEnd;

	__global short* bestQueryPosition;
	__global unsigned char* bestSubjectPosition;
	int bestSegmentScore;
	int nominalScore, count;
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



//Embarrassingly parallel approach is used. One thread is used for
//the hit detection of one sequence
__kernel void search_protein2hitKernel(__global struct PSSMatrixFP *PSSMatrixFP,
										 __global short *matrixBody,
									  	 __global struct sequenceDataFP *sequenceDataFP, 
										 __global unsigned char *sequence,
										 __global struct parameters *parametersFP,
										 __global struct groupFP *wordLookupDFA_groupsFP,
										 __global unsigned char *wordLookupDFAFP,
										 __global uint *blast_numUngappedExtensions,
										 __global uint *blast_numTriggerExtensions,
										 __global uint *blast_numHits,
										 __global uint *hitMatrix_furthestp,
										 __global uint *hitMatrix_offsetp,
										 __global struct ungappedExtension * ungappedExtension_extensionsp,
										 const uint nTotalSequenceNum,
										__global uint* global_numAdditionalTriggerExtensions)
{

int tid = get_global_id(0);

if (tid == 0) {
//printf("GPU ADDRESSES %x %x %x %x %x %x %x %x %x\n\t%x %x %x %x %x %x %x\n", PSSMatrixFP, &PSSMatrixFP[1], matrixBody, sequenceDataFP, sequence, parametersFP, wordLookupDFA_groupsFP, wordLookupDFAFP,
//blast_numUngappedExtensions, blast_numTriggerExtensions, blast_numHits, hitMatrix_furthestp, hitMatrix_offsetp, ungappedExtension_extensionsp, nTotalSequenceNum, global_numAdditionalTriggerExtensions);
//int itr = 0;
//for (; itr < 137072; itr++) {
//	if (hitMatrix_furthestp[itr] != 0) printf("ANGER!! %x %d\n", &hitMatrix_furthestp[itr], hitMatrix_furthestp[itr]);
//}
}

//printf("HELLO from %d!\n", tid);
//	int tid = get_group_id(0) * get_local_size(0) + get_local_id(0);
//	int tid = blockIdx * blockDim + threadIdx;

	__global unsigned char *subject, *sequenceEnd, *address;
	int subjectOffset, count;
	unsigned char currentWord;
	__global unsigned char *currentBlock;
	__global struct groupFP *currentGroupFP;
	__global ushort *wordLookupDFA_AddiPositions;
	uint numOfTriggerExtensions = 0;
	__global ushort *queryOffsets;
	ushort queryOffset;
//	printf("tid: %d\t queryOffset %d\n", tid, queryOffset);
	__global struct ungappedExtension* ungappedExtension_current;
	__global struct ungappedExtension* ungappedExtension_additional;
	int diagonal;
	__global uint *lastHitFP;
	uint ungappedExtension_subjectEndReachedFP;
	__global uint *hitMatrix_Local;
	uint sequenceCount;
	int  distance;


//blast_numHits[tid] = 1000-tid;
//return;

//blast_numTriggerExtensions[tid]=tid;
//ungappedExtension_extensionsp[tid].status = -99;

	hitMatrix_Local = hitMatrix_furthestp +  hitMatrix_offsetp[tid] + PSSMatrixFP->length;
	//ungappedExtension_current = ungappedExtension_extensionsp + tid * UNGAPEXT_PER_THREAD;
	ungappedExtension_extensionsp->start.subjectOffset = 0;
	ungappedExtension_current = ungappedExtension_extensionsp + tid * parametersFP->ungappedExtensionsPerThread;
	ungappedExtension_additional = ungappedExtension_extensionsp + parametersFP->ungappedExtAdditionalStartLoc;
	wordLookupDFA_AddiPositions = (__global unsigned short *)((__global char *)wordLookupDFAFP + 
								  parametersFP->additionalQueryPositionOffset);

	//Set the PSSMatrix body
	PSSMatrixFP->matrix = matrixBody + parametersFP->encoding_numCodes;


int grr = 0;
while (grr < 1){
//printf("brawr %d %d\n", tid, nTotalSequenceNum);
//printf("chargroups[%d] %d %x, groups[%d] %d %x %d\n", grr, ((__global char*)wordLookupDFA_groupsFP)[grr], &((__global char*)wordLookupDFA_groupsFP)[grr], grr, wordLookupDFA_groupsFP[grr].nextWords, &(wordLookupDFA_groupsFP[grr].nextWords), sizeof(int *));
grr++;
}

	sequenceCount = tid;
	while (sequenceCount < nTotalSequenceNum)
	{
		subject = address = sequence + sequenceDataFP[sequenceCount].offset;
	
//printf("addr %d\n", address);
	
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
//				printf("tid %d ", tid);
//				printf("block %d ", currentBlock);
//				printf("groupFP %d ", currentGroupFP);
//				printf("word %d ", currentWord);
//				printf("nextwords %d ", currentGroupFP->nextWords);
//				printf("diff %d ", &wordLookupDFAFP[currentGroupFP->nextWords]-currentBlock);
//				printf("seqEnd %d ", sequenceEnd);
//printf("addr %d\n", address);
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
					queryOffsets = ((__global ushort*)currentBlock) - currentWord;
					
					if (!(*queryOffsets))
					{
						// Go to an outside address for additional positions
						queryOffsets = wordLookupDFA_AddiPositions
									+ ((queryOffsets[1]) * constants_max_short) + queryOffsets[2];
				//					+ (*(queryOffsets + 1) * constants_max_short) + (*(queryOffsets + 2));
					}

					do
					{
						queryOffset = *queryOffsets;
						
//	printf("tid: %d\t queryOffset %d\n", tid, queryOffset);
		
						#ifndef NO_STAGE2
						// Calculate the diagonal this hit is on
						diagonal = subjectOffset - queryOffset;

//printf("diagonal %d %d %d ", diagonal, subjectOffset, queryOffset);
						// If we have not extended past this point on this diagonal
						lastHitFP = hitMatrix_Local + diagonal;
//printf("lastHitFP: %d\n", lastHitFP);

				//		distance = (int)((address - sequence) - *lastHitFP);
						distance = (address - sequence) - *lastHitFP;
//printf("%d ", distance);

if (tid < 167 && tid > 157) {
//printf ("reached tid:%d hit%x diag:%x dist:%x addr:%x seq:%x last:%x lastp:%x parmA:%d parmO:%d\n", tid, hitMatrix_Local, diagonal, distance, address, sequence, *lastHitFP, lastHitFP, parametersFP->parameters_A, parametersFP->parameters_overlap);
}

						if (distance >= parametersFP->parameters_A)
						{	
							//printf("if");
							*lastHitFP = address - sequence; //VUVE
							if (address-sequence > constants_max_short || address-sequence < 0) {
							//	printf("If Violator Assigned from TID: %d!!\n", tid);
							}
						}
						else if (distance >= parametersFP->parameters_overlap)

						{
							//printf("else");
							//Number of extensions for each subject sequence
							blast_numUngappedExtensions[tid]++;

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
									tid,
									global_numAdditionalTriggerExtensions);

							// Update furthest reached value for the diagonal
							*lastHitFP = ungappedExtension_subjectEndReachedFP; //VUVE
							if (ungappedExtension_subjectEndReachedFP > constants_max_short || ungappedExtension_subjectEndReachedFP < 0) {
							//	printf("Else Violator Assigned from TID: %d!!\n", tid);
							}
						}
						#endif
					
						queryOffsets++; blast_numHits[tid]++;
	//		printf("%d %d %d\n", queryOffsets, blast_numUngappedExtensions[tid], numOfTriggerExtensions);
					}
					while((*queryOffsets));
				}
				address++;
			}
		}

		//option=======================================================
		//sequenceCount = atomicAdd(&global_sequenceCount, 1);
		sequenceCount += get_num_groups(0) * get_local_size(0);
//		sequenceCount += gridDim * blockDim;
		//============================================================
	}
blast_numTriggerExtensions[tid] = (uint) numOfTriggerExtensions; mem_fence(CLK_LOCAL_MEM_FENCE);
//printf("tid: %d numTriggerExtensions: %d\n", tid, blast_numTriggerExtensions[tid]);
//	blast_numHits[tid] = tid;

//blast_numTriggerExtensions[tid] = ungappedExtension_extensionsp + tid * parametersFP->ungappedExtensionsPerThread;
//if (tid == 0) blast_numTriggerExtensions[tid] = sequenceCount;
	return;
}
