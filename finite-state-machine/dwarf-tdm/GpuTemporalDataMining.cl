// Implementation copied from glibc
#include "strncpy.cl"
#include "strncmp.cl"
#include "types.h"


//FROM GLOBAL.H
enum
{
	ABSOLUTE,
	RATIO,
	STATIC,
	DYNAMIC,
	MAP_AND_MERGE,
	NAIVE,
	OPTIMAL,
	ONE_CHAR,
	TWO_CHAR
};

#define SDATA( index)      CUT_BANK_CHECKER(sdata, index)
#define FOURBYTEBLOCKINGSIZE 1024

#define MAX_TIMESTAMP_PER_LEVEL 11
#define MAX_LEVEL 20
#define EPSILON 0.0000f
#define IMAGE_MAX_WIDTH 4096
#define IMAGE_POS(length) ((int2)(length % IMAGE_MAX_WIDTH, length / IMAGE_MAX_WIDTH))

int ffs(int x)
{
    int count;
    int var = 1;
    for(count = 1; count <= 32; count++)
    {
        if(x & var)
            return count;
        var = var << 1;
    }
    return 0;
}

void addToList(__local float* values, unsigned int* sizes, unsigned int listIdx, float newValue )
{
	values[listIdx*MAX_TIMESTAMP_PER_LEVEL + sizes[listIdx]] = newValue;
	sizes[listIdx]++;
}

float getFromList(__local float* values, unsigned int listIdx, unsigned int valueIdx )
{
	return values[listIdx*MAX_TIMESTAMP_PER_LEVEL+valueIdx];
}

void clearLists( unsigned int* sizes, const int level )
{
	for ( unsigned int idx = 0; idx < level; idx++ )
		sizes[idx] = 0;
}

////////////////////////////////////////////////////////////////////////////////
//! Simple test kernel for device functionality
//! @param g_idata  input data in global memory
//! @param g_odata  output data in global memory
////////////////////////////////////////////////////////////////////////////////
__kernel void
countCandidates(__global float* episodeSupport, long eventSize, int level, int sType, int numCandidates,
	__read_only image2d_t candidateTex, __read_only image2d_t intervalTex,
	__read_only image2d_t eventTex, __read_only image2d_t timeTex, __local float* timestamps )
{
	unsigned int timestampMatrixSize = level*MAX_TIMESTAMP_PER_LEVEL;
	__local float* timestampMatrix = &timestamps[timestampMatrixSize*get_local_id(0)];
	unsigned int timestampSize[MAX_LEVEL];

	//if ( (get_local_id(0)+get_group_id(0)*get_local_size(0)) >= numCandidates )
	if(get_global_id(0) >= numCandidates)
        return;

	int count = 0;
	int myCandidateIdx = level*(get_global_id(0));//get_local_id(0) + get_group_id(0)*get_local_size(0));
	int myIntervalIdx = 2*(level-1)*(get_global_id(0));//get_local_id(0) + get_group_id(0)*get_local_size(0));

	bool breakOuterLoop;
	// Initialize sizes of timestamp arrays
	clearLists( timestampSize, level );

	for ( long eventIdx = 0; eventIdx < eventSize; eventIdx++ )
	{
		const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
		                          CLK_ADDRESS_CLAMP_TO_EDGE |
								  CLK_FILTER_NEAREST;
		//UBYTE eventSymbol = tex1Dfetch( eventTex, eventIdx );
		UBYTE eventSymbol = read_imageui(eventTex, sampler, IMAGE_POS(eventIdx)).x;
		breakOuterLoop = false;

		// Check other symbols in the episode for matches to the current event
		for ( int symbolIdx = level-1; symbolIdx >= 0 && !breakOuterLoop; symbolIdx-- )
		{
			//if ( eventSymbol == tex1Dfetch( candidateTex, myCandidateIdx + symbolIdx ) )
			if ( eventSymbol ==  read_imageui(candidateTex, sampler, IMAGE_POS(myCandidateIdx + symbolIdx)).x )
			{
				if ( symbolIdx == 0 )
				{
					if ( symbolIdx == level-1 )
						count++;
					else
					{
						//if ( tex1Dfetch(timeTex, eventIdx) - getFromList(timestampMatrix, 0, timestampSize[0]-1) > tex1Dfetch(intervalTex, myIntervalIdx+1) )
						if ( read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x - getFromList(timestampMatrix, 0, timestampSize[0]-1)  >
								read_imagef(intervalTex, sampler, IMAGE_POS(myIntervalIdx+1)).x )
							timestampSize[0] = 0;

						//addToList(timestampMatrix, timestampSize, symbolIdx, tex1Dfetch(timeTex, eventIdx));
						addToList(timestampMatrix, timestampSize, symbolIdx, read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x);
					}
				}
				else
				{
					//float distance = tex1Dfetch(timeTex, eventIdx) - getFromList(timestampMatrix, symbolIdx-1, timestampSize[symbolIdx-1]-1);
					//float lowerBound = tex1Dfetch(intervalTex, myIntervalIdx + (symbolIdx-1)*2+0);
					//float upperBound = tex1Dfetch(intervalTex, myIntervalIdx + (symbolIdx-1)*2+1);
					float distance = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x - getFromList(timestampMatrix, symbolIdx-1, timestampSize[symbolIdx-1]-1);
					float lowerBound = read_imagef(intervalTex, sampler, IMAGE_POS(myIntervalIdx + (symbolIdx-1)*2+0)).x;
					float upperBound = read_imagef(intervalTex, sampler, IMAGE_POS(myIntervalIdx + (symbolIdx-1)*2+1)).x;

					if ( distance > upperBound )
					{
						// Clear list
						timestampSize[symbolIdx-1] = 0;
					}
					else
					{
						// Check previous for acceptable interval
						for ( int prevIdx = timestampSize[symbolIdx-1]-1; prevIdx >= 0; prevIdx-- )
						{
							//distance = tex1Dfetch(timeTex, eventIdx) - getFromList(timestampMatrix, symbolIdx-1, prevIdx);
							distance = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x - getFromList(timestampMatrix, symbolIdx-1, prevIdx);

							if ( (distance >  lowerBound  - EPSILON &&
								  distance <= upperBound  + EPSILON )||
								  symbolIdx == 0)
							{
								if ( symbolIdx == level-1 )
								{
									// The final symbol has been found, clear all lists
									count++;
									clearLists( timestampSize, level );
									breakOuterLoop = true;
								}
								else
								{
									//if ( tex1Dfetch(timeTex, eventIdx) - getFromList(timestampMatrix, symbolIdx, timestampSize[symbolIdx]-1) > tex1Dfetch(intervalTex, myIntervalIdx + 2*(symbolIdx)+1) )
									if ( read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x - getFromList(timestampMatrix, symbolIdx, timestampSize[symbolIdx]-1) > read_imagef(intervalTex, sampler, IMAGE_POS(myIntervalIdx + 2*(symbolIdx)+1)).x )
										timestampSize[symbolIdx] = 0;

									addToList(timestampMatrix, timestampSize, symbolIdx, read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x);
								}
								break;
							}
						}
					}
				}
			}
		}
	}


	episodeSupport[get_local_id(0) + get_group_id(0)*get_local_size(0)] = sType == RATIO ? ((float)count / (float)eventSize) : (float)count;
}

void resetToZero(__local float* timestamps, int level )
{
	for( int idx = 0; idx < level; idx++ )
		timestamps[idx] = -1.0f;
}

__kernel void
countCandidatesStatic(__global float* episodeSupport, long eventSize, int level, int sType, int numCandidates,
	__read_only image2d_t candidateTex, __read_only image2d_t intervalTex,
	__read_only image2d_t eventTex, __read_only image2d_t timeTex,__local float* timestamps)
{
	__local float* myTimestamps = &timestamps[get_local_id(0)*level];

	int count = 0;
	int myCandidateIdx = level*(get_local_id(0) + get_group_id(0)*get_local_size(0));
	int myIntervalIdx = 2*(level-1)*(get_local_id(0) + get_group_id(0)*get_local_size(0));

	if ( (get_local_id(0)+get_group_id(0)*get_local_size(0)) >= numCandidates )
		return;

	// Initialize sizes of timestamp arrays
	resetToZero( myTimestamps, level );

	for ( long eventIdx = 0; eventIdx < eventSize; eventIdx++ )
	{
		const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
		                          CLK_ADDRESS_NONE |
								  CLK_FILTER_NEAREST;
		//UBYTE eventSymbol = tex1Dfetch( eventTex, eventIdx );
		UBYTE eventSymbol = read_imageui( eventTex, sampler, IMAGE_POS(eventIdx)).x;

		// Check other symbols in the episode for matches to the current event
		for ( int symbolIdx = level-1; symbolIdx >= 0; symbolIdx-- )
		{
			//if ( eventSymbol == tex1Dfetch( candidateTex, myCandidateIdx + symbolIdx ) )
			if ( eventSymbol == read_imageui( candidateTex, sampler, IMAGE_POS(myCandidateIdx + symbolIdx)).x )
			{
				if ( symbolIdx == 0 )
				{
					if ( symbolIdx == level-1 )
						count++;
					else
					{
						//myTimestamps[symbolIdx] = tex1Dfetch(timeTex, eventIdx);
						myTimestamps[symbolIdx] = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x;					
					}
				}
				else
				{
					//float distance = tex1Dfetch(timeTex, eventIdx) - myTimestamps[symbolIdx-1];
					//float upperBound = tex1Dfetch(intervalTex, myIntervalIdx + (symbolIdx-1)*2+1);
					float distance = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x - myTimestamps[symbolIdx-1];
					float upperBound = read_imagef(intervalTex, sampler, IMAGE_POS(myIntervalIdx + (symbolIdx-1)*2+1)).x;

					if (  distance <= upperBound + EPSILON ||
						  symbolIdx == 0)
					{
						if ( symbolIdx == level-1 )
						{
							// The final symbol has been found, clear all lists
							count++;
							resetToZero(myTimestamps, level);
							break;
						}
						else
						{
							//myTimestamps[symbolIdx] = tex1Dfetch(timeTex, eventIdx);
							myTimestamps[symbolIdx] = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x;
						}
					}
				}
			}
		}
	}

	episodeSupport[get_local_id(0) + get_group_id(0)*get_local_size(0)] = sType == RATIO ? ((float)count / (float)eventSize) : (float)count;
}

__kernel
void countCandidatesMapMerge(__global float* episodeSupport, long eventSize, int level, int sType, int numSections, int eventsPerSection, int numCandidates,
	__read_only image2d_t candidateTex, __read_only image2d_t intervalTex,
	__read_only image2d_t eventTex, __read_only image2d_t timeTex, __local float* timestamps)
{
	unsigned int timestampMatrixSize = (level-1)*MAX_TIMESTAMP_PER_LEVEL;
	__local float* timestampMatrix = &timestamps[timestampMatrixSize*(get_local_id(1)*get_local_size(0)+get_local_id(0))];
	unsigned int timestampSize[MAX_TIMESTAMP_PER_LEVEL];
	__local float* records = &timestamps[level*numSections*timestampMatrixSize];

    const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
							  CLK_ADDRESS_NONE |
							  CLK_FILTER_NEAREST;
								  
	int count = 0;
	int tid = get_local_id(1)*get_local_size(0) + get_local_id(0);
	int bid = get_group_id(1)*get_num_groups(0) + get_group_id(0);

	if ( bid > numCandidates )
		return;

	int myCandidateIdx = level*bid;
	int myIntervalIdx = 2*(level-1)*bid;
	float episodeLength = 0.0f;
	float partialEpisodeLength = 0.0f;

	int firstOccurrence = -1, lastOccurrence = -1;

	// Calculate first and last element in section
	int sectionStart = get_local_id(0)*eventsPerSection;
	int sectionEnd = ((get_local_id(0)+1)*eventsPerSection >= eventSize) ? eventSize-1 : (get_local_id(0)+1)*eventsPerSection-1;

	// Find Max Length of Episode
	for ( unsigned int idx = 0; idx < level-1; idx++ )
	{
		// Resetting timestamp size since we already have a loop going
		timestampSize[idx] = 0;

		episodeLength += read_imagef( intervalTex, sampler, IMAGE_POS( myIntervalIdx + 2*idx+1) ).x;

		// Find Length for my state machine
		if ( idx < get_local_id(1))
			partialEpisodeLength = episodeLength;
	}

	float startTime = read_imagef( timeTex, sampler, IMAGE_POS( sectionStart)).x - partialEpisodeLength;
	float endTime = read_imagef( timeTex, sampler, IMAGE_POS(sectionEnd) ).x + episodeLength;


	// Determine start index
	int startIdx = sectionStart;
	while( startIdx > 0 && read_imagef( timeTex, sampler, IMAGE_POS( startIdx) ).x > startTime )
		startIdx--;

	bool breakOuterLoop;

	// Initialize sizes of timestamp arrays
	clearLists( timestampSize, level-1 );

	for ( long eventIdx = startIdx; eventIdx < eventSize && read_imagef( timeTex, sampler, IMAGE_POS(eventIdx) ).x <= endTime; eventIdx++ )
	{
		UBYTE eventSymbol = read_imageui( eventTex, sampler, IMAGE_POS(eventIdx)).x;
		breakOuterLoop = false;

		// Check other symbols in the episode for matches to the current event
		for ( int symbolIdx = level-1; symbolIdx >= 0 && !breakOuterLoop; symbolIdx-- )
		{
			UBYTE searchSymbol = read_imageui( candidateTex, sampler, IMAGE_POS( myCandidateIdx + symbolIdx) ).x;

//			if ( searchSymbol == 'Q' )
//				searchSymbol = 'Q';

			//if ( eventSymbol == tex1Dfetch( candidateTex, myCandidateIdx + symbolIdx ) )
			if ( eventSymbol == searchSymbol )
			{
				if ( symbolIdx == 0 )
				{
					if ( symbolIdx == level-1 )
					{
						if ( eventIdx >= sectionStart - EPSILON && eventIdx <= sectionEnd + EPSILON )
							count++;

						// Record first and last occurrences
						if ( firstOccurrence == -1 )
							firstOccurrence = eventIdx;

						lastOccurrence = eventIdx;
					}
					else
					{
						if ( timestampSize[0] != 0 && read_imagef(timeTex, sampler, IMAGE_POS( eventIdx)).x - getFromList(timestampMatrix, 0, timestampSize[0]-1) > read_imagef(intervalTex, sampler, IMAGE_POS( myIntervalIdx+1)).x )
							timestampSize[0] = 0;

						addToList(timestampMatrix, timestampSize, symbolIdx, read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x);
					}
				}
				else
				{
					if ( timestampSize[symbolIdx-1] == 0 )
						continue;
					float distance = read_imagef(timeTex, sampler, IMAGE_POS( eventIdx)).x - getFromList(timestampMatrix, symbolIdx-1, timestampSize[symbolIdx-1]-1);
					float lowerBound = read_imagef(intervalTex, sampler, IMAGE_POS( myIntervalIdx + (symbolIdx-1)*2+0)).x;
					float upperBound = read_imagef(intervalTex, sampler, IMAGE_POS( myIntervalIdx + (symbolIdx-1)*2+1)).x;

					if ( distance > upperBound )
					{
						// Clear list
						timestampSize[symbolIdx-1] = 0;
					}
					else
					{
						// Check previous for acceptable interval
						for ( int prevIdx = timestampSize[symbolIdx-1]-1; prevIdx >= 0 && timestampSize[symbolIdx-1] != 0; prevIdx-- )
						{
							//distance = tex1Dfetch(timeTex, eventIdx) - getFromList(timestampMatrix, symbolIdx-1, prevIdx);
							distance = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x - getFromList(timestampMatrix, symbolIdx-1, prevIdx);

							if ( (distance >  lowerBound  - EPSILON &&
								  distance <= upperBound  + EPSILON )||
								  symbolIdx == 0)
							{
								if ( symbolIdx == level-1 )
								{
									// The final symbol has been found, increment and clear all lists
									if ( eventIdx >= sectionStart - EPSILON && eventIdx <= sectionEnd + EPSILON )
										count++;

									breakOuterLoop = true;
									clearLists( timestampSize, level-1 );

									// Record first and last occurrences
									if ( firstOccurrence == -1 )
										firstOccurrence = eventIdx;

									lastOccurrence = eventIdx;
								}
								else
								{
									//float ts = tex1Dfetch(timeTex, eventIdx);
									//float last = getFromList(timestampMatrix, symbolIdx, timestampSize[symbolIdx]-1);
									//float upper = tex1Dfetch(intervalTex, myIntervalIdx + 2*(symbolIdx+1)+1);

									//if ( timestampSize[symbolIdx] != 0 && tex1Dfetch(timeTex, eventIdx) - getFromList(timestampMatrix, symbolIdx, timestampSize[symbolIdx]-1) > tex1Dfetch(intervalTex, myIntervalIdx + 2*(symbolIdx)+1) )
									if ( timestampSize[symbolIdx] != 0 && read_imagef(timeTex, sampler, IMAGE_POS( eventIdx)).x - getFromList(timestampMatrix, symbolIdx, timestampSize[symbolIdx]-1) > read_imagef(intervalTex, sampler, IMAGE_POS( myIntervalIdx + 2*(symbolIdx)+1) ).x)
									//if ( timestampSize[symbolIdx] != 0 && ts - last > upper )
										timestampSize[symbolIdx] = 0;

									//addToList(timestampMatrix, timestampSize, symbolIdx, tex1Dfetch(timeTex, eventIdx));
									addToList(timestampMatrix, timestampSize, symbolIdx, read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x);
								}
								break;
							}
						}
					}
				}
			}
		}
	}

	// Write first, last occurrence and count in shared memory
	records[3*level*get_local_id(0)+3*get_local_id(1)+0] = firstOccurrence;
	records[3*level*get_local_id(0)+3*get_local_id(1)+1] = lastOccurrence;
	records[3*level*get_local_id(0)+3*get_local_id(1)+2] = count;

	barrier(CLK_LOCAL_MEM_FENCE|CLK_GLOBAL_MEM_FENCE);

	//start reduction here
	int mergeScope = numSections;
	int mergeSpace = 1;
	int bestMatch;
	//#pragma unroll
	for ( int it = 0; it < ffs(numSections) - 1; it++ )
	{
		mergeScope = mergeScope >> 1;
		if ( tid < mergeScope )
		{
			int first = 2*mergeSpace*tid;
			int second = first + mergeSpace;

			// Find Matching
			for ( int leftIdx = 0; leftIdx < level; leftIdx++ )
			{
				bestMatch = -1;
				for ( int rightIdx = 0; rightIdx < level; rightIdx++ )
				{
					// Match left->last to right->first
					if ( records[3*level*first + 3*leftIdx + 1] == records[3*level*second + 3*rightIdx + 0] )
					{
						bestMatch = rightIdx;
					}
					else if ( (records[3*level*first + 3*leftIdx + 1] < records[3*level*second + 3*rightIdx + 0]
							|| records[3*level*second + 3*rightIdx + 0] == -1
							|| records[3*level*first + 3*leftIdx + 1] == -1)
								&& bestMatch == -1 )
					{
						bestMatch = rightIdx;
					}
				}

				// Set last occurrence to best match's last occurrence, add counts
				records[3*level*first + 3*leftIdx + 1] = records[3*level*second + 3*bestMatch + 1];
				records[3*level*first + 3*leftIdx + 2] += records[3*level*second + 3*bestMatch + 2];
			}

		}
		barrier(CLK_LOCAL_MEM_FENCE|CLK_GLOBAL_MEM_FENCE);

		mergeSpace = mergeSpace << 1;
	}



























	//int width = 2;
	//int bestMatch;
	//bool mergeComplete = false;
	//for ( int i = (blockDim.y*get_local_size(0))/2; !mergeComplete; i/=2 )
	//{
	//	if ( i == 0 )
	//		mergeComplete = true;

	//	if ( tid <= i )
	//	{
	//		// Find Matching
	//		for ( int leftIdx = 0; leftIdx < level; leftIdx++ )
	//		{
	//			bestMatch = -1;
	//			for ( int rightIdx = 0; rightIdx < level; rightIdx++ )
	//			{
	//				// Match left->last to right->first
	//				if ( records[3*level*tid*width + 3*leftIdx + 1] == records[3*level*(tid*width + width/2) + 3*rightIdx + 0] )
	//				{
	//					bestMatch = rightIdx;
	//				}
	//				else if ( records[3*level*tid*width + 3*leftIdx + 1] < records[3*level*(tid*width + width/2) + 3*rightIdx + 0]
	//							&& bestMatch == -1 )
	//				{
	//					bestMatch = rightIdx;
	//				}
	//			}

	//			// Set last occurrence to best match's last occurrence, add counts
	//			records[3*level*tid*width + 3*leftIdx + 1] = records[3*level*(tid*width + width/2) + 3*bestMatch + 1];
	//			records[3*level*tid*width + 3*leftIdx + 2] += records[3*level*(tid*width + width/2) + 3*bestMatch + 2];
	//		}
	//	}
	//
	//	width *= 2;
	//	__syncthreads();
	//}

	//__syncthreads();

	// Get max value
	if ( tid == 0 )
	{
		count = records[2];
		for ( unsigned int idx = 1; idx < level; idx++ )
		{
			if ( count < records[3*idx+2] )
				count = records[3*idx+2];
		}

		episodeSupport[bid] = sType == RATIO ? ((float)count / (float)eventSize) : (float)count;
	}
}




__kernel
void countCandidatesMapMergeStatic(__global float* episodeSupport, long eventSize, int level, int sType, int numSections, int eventsPerSection, int numCandidates,
	__read_only image2d_t candidateTex, __read_only image2d_t intervalTex,
	__read_only image2d_t eventTex, __read_only image2d_t timeTex, __local float* timestamps )
{

	const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
							  CLK_ADDRESS_NONE |
							  CLK_FILTER_NEAREST;
	int count = 0;
	int tid = get_local_id(1)*get_local_size(0) + get_local_id(0);
	int bid = get_group_id(1)*get_num_groups(0) + get_group_id(0);

	if ( bid > numCandidates )
		return;

	int myCandidateIdx = level*bid;
	int myIntervalIdx = 2*(level-1)*bid;
	float episodeLength = 0.0f;
	float partialEpisodeLength = 0.0f;

	__local float* myTimestamps = &timestamps[tid*level];
	__local float* records = &timestamps[level*level*numSections];

	int firstOccurrence = -1, lastOccurrence = -1;

	// Calculate first and last element in section
	int sectionStart = get_local_id(0)*eventsPerSection;
	int sectionEnd = ((get_local_id(0)+1)*eventsPerSection >= eventSize) ? eventSize-1 : (get_local_id(0)+1)*eventsPerSection-1;

	// Find Max Length of Episode
	for ( unsigned int idx = 0; idx < level-1; idx++ )
	{
		//episodeLength += tex1Dfetch( intervalTex, myIntervalIdx + 2*idx+1 );
		episodeLength += read_imagef( intervalTex, sampler, IMAGE_POS( myIntervalIdx + 2*idx+1) ).x;

		// Find Length for my state machine
		if ( idx < get_local_id(1) )
			partialEpisodeLength = episodeLength;
	}

	float startTime = read_imagef( timeTex, sampler, IMAGE_POS( sectionStart )).x - partialEpisodeLength;
	float endTime = read_imagef( timeTex, sampler, IMAGE_POS( sectionEnd )).x + episodeLength;


	// Determine start index
	int startIdx = sectionStart;
	//while( startIdx > 0 && tex1Dfetch( timeTex, startIdx ) >= startTime )
	while( startIdx > 0 && read_imagef(timeTex, sampler, IMAGE_POS( startIdx) ).x >= startTime )
		startIdx--;

	// Initialize sizes of timestamp arrays
	resetToZero( myTimestamps, level );

	//for ( long eventIdx = startIdx; eventIdx < eventSize && tex1Dfetch( timeTex, eventIdx ) <= endTime; eventIdx++ )
	for ( long eventIdx = startIdx; eventIdx < eventSize && read_imagef( timeTex, sampler, IMAGE_POS( eventIdx)).x  <= endTime; eventIdx++ )
	{
		//UBYTE eventSymbol = tex1Dfetch( eventTex, eventIdx );
		UBYTE eventSymbol = read_imageui( eventTex, sampler, IMAGE_POS(eventIdx) ).x;

		// Check other symbols in the episode for matches to the current event
		for ( int symbolIdx = level-1; symbolIdx >= 0; symbolIdx-- )
		{
			//if ( eventSymbol == tex1Dfetch( candidateTex, myCandidateIdx + symbolIdx ) )
			if ( eventSymbol == read_imageui( candidateTex, sampler, IMAGE_POS( myCandidateIdx + symbolIdx)).x )
			{
				if ( symbolIdx == 0 )
				{
					if ( symbolIdx == level-1 )
					{
						if ( eventIdx >= sectionStart && eventIdx <= sectionEnd )
							count++;

						// Record first and last occurrences
						if ( firstOccurrence == -1 )
							firstOccurrence = eventIdx;

						lastOccurrence = eventIdx;
					}
					else
					{
						//myTimestamps[symbolIdx] = tex1Dfetch(timeTex, eventIdx);
						myTimestamps[symbolIdx] = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x;					
					}
				}
				else
				{
					//float distance = tex1Dfetch(timeTex, eventIdx) - myTimestamps[symbolIdx-1];
					//float upperBound = tex1Dfetch(intervalTex, myIntervalIdx + (symbolIdx-1)*2+1);
					float distance = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x - myTimestamps[symbolIdx-1];
					float upperBound = read_imagef(intervalTex, sampler, IMAGE_POS( myIntervalIdx + (symbolIdx-1)*2+1)).x;

					if (  distance <= upperBound + EPSILON ||
						  symbolIdx == 0)
					{
						if ( symbolIdx == level-1 )
						{
							// The final symbol has been found, increment and clear all lists
							if ( eventIdx >= sectionStart && eventIdx <= sectionEnd )
								count++;

							// Record first and last occurrences
							if ( firstOccurrence == -1 )
								firstOccurrence = eventIdx;

							lastOccurrence = eventIdx;

							// The final symbol has been found, clear all lists
							resetToZero(myTimestamps, level);
							break;
						}
						else
						{
							//myTimestamps[symbolIdx] = tex1Dfetch(timeTex, eventIdx);
							myTimestamps[symbolIdx] = read_imagef(timeTex, sampler, IMAGE_POS(eventIdx)).x;
						}
					}
				}
			}
		}
	}

	// Write first, last occurrence and count in shared memory
	records[3*level*get_local_id(0)+3*get_local_id(1)+0] = firstOccurrence;
	records[3*level*get_local_id(0)+3*get_local_id(1)+1] = lastOccurrence;
	records[3*level*get_local_id(0)+3*get_local_id(1)+2] = count;

	barrier(CLK_LOCAL_MEM_FENCE|CLK_GLOBAL_MEM_FENCE);


	//start reduction here
	int mergeScope = numSections;
	int mergeSpace = 1;
	int bestMatch;
	//#pragma unroll
	for ( int it = 0; it < ffs(numSections) - 1; it++ )
	{
		mergeScope = mergeScope >> 1;
		if ( tid < mergeScope )
		{
			int first = 2*mergeSpace*tid;
			int second = first + mergeSpace;

			// Find Matching
			for ( int leftIdx = 0; leftIdx < level; leftIdx++ )
			{
				bestMatch = -1;
				for ( int rightIdx = 0; rightIdx < level; rightIdx++ )
				{
					// Match left->last to right->first
					if ( records[3*level*first + 3*leftIdx + 1] == records[3*level*second + 3*rightIdx + 0] )
					{
						bestMatch = rightIdx;
					}
					else if ( (records[3*level*first + 3*leftIdx + 1] < records[3*level*second + 3*rightIdx + 0]
							|| records[3*level*second + 3*rightIdx + 0] == -1
							|| records[3*level*first + 3*leftIdx + 1] == -1)
								&& bestMatch == -1 )
					{
						bestMatch = rightIdx;
					}
				}

				// Set last occurrence to best match's last occurrence, add counts
				records[3*level*first + 3*leftIdx + 1] = records[3*level*second + 3*bestMatch + 1];
				records[3*level*first + 3*leftIdx + 2] += records[3*level*second + 3*bestMatch + 2];
			}

		}
		barrier(CLK_LOCAL_MEM_FENCE|CLK_GLOBAL_MEM_FENCE);

		mergeSpace = mergeSpace << 1;
	}

	// Get max value
	if ( tid == 0 )
	{
		count = records[2];
		for ( unsigned int idx = 1; idx < level; idx++ )
		{
			if ( count < records[3*idx+2] )
				count = records[3*idx+2];
		}

		episodeSupport[bid] = sType == RATIO ? ((float)count / (float)eventSize) : (float)count;
	}
}







//__device__ void
//initEpisodeCandidates(UBYTE* episodeCandidates)
//{
//	if ( get_local_id(0) < uniqueEvents )
//	{
//		episodeCandidates[get_local_id(0)] = get_local_id(0) + 'A';
//	}
//}

// Convert to GPU code
//__device__ int
//getNextValidCandidateGPU(int prefixLength, int currentIdx, int nextIdx, int numCandidates,
//						UBYTE* episodeCandidates, float* episodeSupport)
//{
//	nextIdx++;
//	for ( int idx = nextIdx; idx < numCandidates; idx++)
//	{
//		if ( cuda_strncmp( (char*)&episodeCandidates[currentIdx*(prefixLength+1)],
//					  (char*)&episodeCandidates[idx*(prefixLength+1)],
//					  prefixLength ) == 0 && episodeSupport[idx] > 0.0f )
//		{
//			return idx;
//		}
//	}
//
//	return -1;
//}
//

// Converts "triangle" matrix to two array positions - the base position, and compare position
void
triangleToArray( int triangle, int* base, int* compare, int numCandidates )
{
	*base = 0;
	int Temp = triangle + 1;
	while (Temp > 0)
	{
   		Temp = Temp - (numCandidates - *base - 1);
   		(*base)++;
	}

	Temp += numCandidates - *base - 1;
	*compare = Temp + (*base);
	(*base)--;
}

bool
compareCandidates( __global UBYTE* episodeCandidates, int level, int base, int compare )
{
	int idx;
	for ( idx = 0; idx < level-2; idx++ )
	{
		if ( episodeCandidates[base*(level-1) + idx] != episodeCandidates[compare*(level-1) + idx] )
		{
			return false;
		}
	}

	return true;
}


// INPUT: stuff
// OUTPUT: Triangle array filled with matching pairs
__kernel void
analyzeSupport( float support,  int level, int numPairs,
				__global UBYTE* episodeCandidates,
				__global float* episodeSupport, __global bool* episodePairs,
				int numCandidates )
{
	int triangleIndex = get_group_id(0) * get_local_size(0) + get_local_id(0);

	if ( triangleIndex >= numPairs )
		return;

	int base, compare;	// base is the original element, compare is the element to compare to
	triangleToArray( triangleIndex, &base, &compare, numCandidates );

	if ( episodeSupport[base] < support || episodeSupport[compare] < support )
	{
		episodePairs[triangleIndex] = false;
		return;
	}

	episodePairs[triangleIndex] = compareCandidates( episodeCandidates, level, base, compare );
}

//Generate episodes from level - 1 to level
__kernel void
generateEpisodeCandidatesKernel( int level, int numCandidates, int numCandidatesBuffer, __global UBYTE* episodeCandidates, __global int* episodeIndices, __global UBYTE* episodeCandidatesBuffer )
{
	int episodeIndex = get_group_id(0) * get_local_size(0) + get_local_id(0);

	if ( episodeIndex >= numCandidatesBuffer )
		return;

	int first, second;
	triangleToArray( episodeIndices[episodeIndex], &first, &second, numCandidates );

	int idx;
	int bufferBaseIndex = 2*level*episodeIndex;
	for ( idx = 0; idx < level-2; idx++ )
	{
		episodeCandidatesBuffer[bufferBaseIndex + idx] = episodeCandidates[first*(level-1) + idx];
		episodeCandidatesBuffer[bufferBaseIndex + level + idx] = episodeCandidates[first*(level-1) + idx];
	}
	episodeCandidatesBuffer[bufferBaseIndex + level - 2] = episodeCandidates[first*(level-1) + level - 2];
	episodeCandidatesBuffer[bufferBaseIndex + level - 1] = episodeCandidates[second*(level-1) + level - 2];

	episodeCandidatesBuffer[bufferBaseIndex + 2 * level - 2] = episodeCandidates[second*(level-1) + level - 2];
	episodeCandidatesBuffer[bufferBaseIndex + 2 * level - 1] = episodeCandidates[first*(level-1) + level - 2];

}
