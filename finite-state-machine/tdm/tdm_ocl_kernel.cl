// countCandidates
//
// targetEvent - Event that we're looking for, this should be pushed into an array if there are more than 1 episode candidate
// startRecords - start/end times of episode candidate evidence from the previous level
// foundRecords - output, store start/end times for evidence for this level
//
__kernel void writeCandidates( unsigned char targetEvent, float minInterval, float maxInterval, 
		__global uint2* startRecords, __global unsigned int* startOffset, __global uint2* foundRecords, 
		uint recordCount, __global unsigned char* d_events, __global float* d_times ) {

	uint tid = get_local_size(0)*get_group_id(0)+get_local_id(0);
	uint2 ep;	// .x is start, .y is end
	uint myIdx, offset;
	float startTime, curTime;

	unsigned int listSize = 0;

	if ( tid < recordCount ) {
		// Retrieve endIndex from startRecords and set to startIdx
		ep = startRecords[tid];
		offset = startOffset[tid];
		myIdx = ep.x;

		// Retrieve interval for this event
		curTime = startTime = d_times[myIdx];

		// Advance through events until minimum interval has been reached
		while( startTime - curTime < minInterval && myIdx > 0 )	{
			myIdx--;
			curTime = d_times[myIdx];
		}

		// While the current event has not exceeded the maximum interval
		while( startTime - curTime <= maxInterval && myIdx > 0) {
			if ( d_events[myIdx] == targetEvent ) {
				foundRecords[offset + listSize].x=myIdx;
				foundRecords[offset + listSize].y=ep.y;
				listSize++;
			}
			myIdx--;
			curTime = d_times[myIdx];
		}
	}
}



// countCandidates
//
// targetEvent - Event that we're looking for, this should be pushed into an array if there are more than 1 episode candidate
// startRecords - start/end times of episode candidate evidence from the previous level
// foundRecords - output, store start/end times for evidence for this level
//
__kernel
void countCandidates( unsigned char targetEvent, float minInterval, float maxInterval, 
		__global uint2* startRecords, __global uint* foundRecordsCount, uint recordCount, __global unsigned char* d_events, __global float* d_times ) {

	uint tid = get_local_size(0)*get_group_id(0)+get_local_id(0);
	uint2 ep;	// .x is start, .y is end
	uint myIdx;
	float startTime, curTime;

	unsigned int listSize = 0;

	if ( tid < recordCount ) {
		// Retrieve endIndex from startRecords and set to startIdx
		ep = startRecords[tid];
		myIdx = ep.x;

		// Retrieve interval for this event
		curTime = startTime = d_times[myIdx];

		// Advance through events until minimum interval has been reached
		while( startTime - curTime < minInterval && myIdx > 0 )	{
			myIdx--;
			curTime = d_times[myIdx];
		}

		// While the current event has not exceeded the maximum interval
		while( startTime - curTime <= maxInterval && myIdx > 0) {
			if ( d_events[myIdx] == targetEvent ) {
				listSize++;
			}//if
			myIdx--;
			curTime = d_times[myIdx];
		}//while
		foundRecordsCount[tid] = listSize;
	}//if
}
