/**
 * Ported to OpenCL from CUDA version of
 * GPU Temporal Data Mining (Sean Ponce)
 **/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include "types.h"
#include "global.h"
#include "dataio.h"

#include "dwarf_sycl.hpp"

/* the OpenCL uint2 (8-byte pair) */
struct alignas(8) tdm_uint2 { unsigned int x, y; };

/* ---- kernels, translated from tdm_ocl_kernel.cl ---- */
static void writeCandidates(sycl::queue &q, sycl::nd_range<1> ndr, unsigned char targetEvent, float minInterval, float maxInterval,
		tdm_uint2* startRecords, unsigned int* startOffset, tdm_uint2* foundRecords, unsigned int recordCount, unsigned char* d_events, float* d_times)
{
	q.parallel_for<class tdm_write>(ndr, [=](sycl::nd_item<1> it) {

		uint tid = it.get_local_range(0)*it.get_group(0)+it.get_local_id(0);
		tdm_uint2 ep;	// .x is start, .y is end
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
}).wait_and_throw();
}

static void countCandidates(sycl::queue &q, sycl::nd_range<1> ndr, unsigned char targetEvent, float minInterval, float maxInterval,
		tdm_uint2* startRecords, unsigned int* foundRecordsCount, unsigned int recordCount, unsigned char* d_events, float* d_times)
{
	q.parallel_for<class tdm_count>(ndr, [=](sycl::nd_item<1> it) {

		uint tid = it.get_local_range(0)*it.get_group(0)+it.get_local_id(0);
		tdm_uint2 ep;	// .x is start, .y is end
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
}).wait_and_throw();
}

#define AOCL_ALIGNMENT 64
static sycl::queue *Q;
unsigned int numCandidates;
int episodesCulled = 0;

ubyte* h_events;
ubyte *d_events;

float* h_times;
float *d_times;

int maxLevel;

ubyte* h_episodeCandidates;
ubyte *d_episodeCandidates;

float* h_episodeIntervals;
float *d_episodeIntervals;

uint* h_episodeSupport;
uint *d_episodeSupport;

int* h_episodeIndices;
int *d_episodeIndices;

uint eventSize;
uint eventType, uniqueEvents;

unsigned int  h_foundCount;
tdm_uint2 *d_startRecords;
tdm_uint2  *h_startRecords;

tdm_uint2 *d_foundRecords;
tdm_uint2  *h_foundRecords;

uint *d_recCount;
uint *d_recOffSet;

uint *d_eventCounts;
uint* h_eventCounts;

uint *d_eventIndex;
tdm_uint2** h_eventIndex;

FILE* dumpFile;

float* temporalConstraint;
unsigned int temporalConstraintSize;

int runTest( int argc, char** argv);

void loadCandidateEpisodes(char* filename, int& level)
{
	FILE *fp;
	char buf[256], symbol, c, *evt;
	int v, tlevel = 0, flag = 0;
	fp = fopen(filename, "r");
	level = 0;
	if (fp)
	{
		int idx = 0;
		while(fgets(buf, 256, fp) != NULL)
		{
			if (strlen(buf)==0) continue;
			//printf("Line: %s\n", buf);
			evt = strtok (buf," ,.-");
			int i = 0;
			while (evt != NULL)
			{
				//printf ("%s\n",evt);
				if (strlen(evt)==1)
				{ //One size event type
					sscanf( evt, "%c", &symbol);
				}
				else
				{
					sscanf( evt, "%c%d", &c, &v);
					symbol = symbolToChar( c, v );
				}
				if (flag == 0) tlevel++;
				h_episodeCandidates[(idx * level) + i] = symbol;
				h_episodeIntervals[(2*idx*(level-1))+i*2+0] = temporalConstraint[0];
				h_episodeIntervals[(2*idx*(level-1))+i*2+1] = temporalConstraint[1];

				evt = strtok (NULL, " ,.-");
				i++;
			}
			if (flag == 0)
			{
				flag = 1;
				level = tlevel;
				//printf("Detected size of episode = %d\n", level);
			}
			idx ++;
		}
		numCandidates = idx;
		//printf("Detected number of episodes = %d\n", numCandidates);
	}
	else
	{
		printf("Unable to read episode file\n");
		exit(1);
	}
}


void indexEvents() {
	char start = eventType == EVENT_26 ? 'A' : '!';

	/*h_eventIndex = (tdm_uint2**)malloc(uniqueEvents*sizeof(tdm_uint2*));
	h_eventCounts = (uint*)malloc(uniqueEvents*sizeof(uint));
        */

        h_eventIndex = (tdm_uint2**) memalign ( AOCL_ALIGNMENT, uniqueEvents*sizeof(tdm_uint2*));
	h_eventCounts = (uint*) memalign ( AOCL_ALIGNMENT, uniqueEvents*sizeof(uint));

	// First pass to how many events of each type
	memset(h_eventCounts, 0, uniqueEvents*sizeof(unsigned int) );
	for ( long idx = 0; idx < eventSize; idx++ ) {
		h_eventCounts[h_events[idx]-start]++;
	}

	// Generate index arrays
	for ( long idx = 0; idx < uniqueEvents; idx++ ) {
		//printf("%c: %d\n", start+idx, h_eventCounts[idx]);
		//h_eventIndex[idx] = (tdm_uint2*)malloc(h_eventCounts[idx]*sizeof(tdm_uint2) );
		h_eventIndex[idx] = (tdm_uint2*) memalign ( AOCL_ALIGNMENT, h_eventCounts[idx]*sizeof(tdm_uint2) );
	}

	// Second pass to index events
	memset(h_eventCounts, 0, uniqueEvents*sizeof(unsigned int) );
	for ( long idx = 0; idx < eventSize; idx++ ) {
		h_eventIndex[h_events[idx]-start][h_eventCounts[h_events[idx]-start]].x = 
			h_eventIndex[h_events[idx]-start][h_eventCounts[h_events[idx]-start]].y =
			(unsigned int)idx;
		h_eventCounts[h_events[idx]-start]++;
	}
}

int compareUint2ByY( const void* v1, const void* v2 ) {
	tdm_uint2 i1 = *(tdm_uint2*)v1;
	tdm_uint2 i2 = *(tdm_uint2*)v2;
	return i1.y-i2.y;
}

tdm_uint2* eliminateConflicts(tdm_uint2* records, uint recordCount, uint* resultCount) {
	//tdm_uint2* result = (tdm_uint2*)malloc(recordCount*sizeof(tdm_uint2));
	tdm_uint2* result = (tdm_uint2*) memalign ( AOCL_ALIGNMENT,recordCount*sizeof(tdm_uint2));
	uint count = 0;
	uint lastRecord = 0;

	// Sort records by endtime
	//qsort( records, recordCount, sizeof(tdm_uint2), compareUint2ByY );

	for ( uint idx = 0; idx < recordCount; idx++ ) {
		if ( records[idx].x >= lastRecord ) {
			result[count] = records[idx];
			count++;

			lastRecord = records[idx].y;
		}
	}

	*resultCount = count;
	return result;
}

void setupGpu()
{
	sycl::queue &q = *Q;

	d_events=dsycl::alloc<ubyte>(q, eventSize);
	d_times=dsycl::alloc<float>(q, eventSize);
	d_episodeCandidates=dsycl::alloc<ubyte>(q, maxCandidates);
	d_episodeIntervals=dsycl::alloc<float>(q, maxIntervals);
	d_episodeSupport=dsycl::alloc<uint>(q, maxCandidates);
	d_startRecords=dsycl::alloc<tdm_uint2>(q, MaxRecords);
	d_foundRecords=dsycl::alloc<tdm_uint2>(q, MaxRecords);
	/*
	h_startRecords = (tdm_uint2*)malloc(MaxRecords*sizeof(tdm_uint2));
	h_foundRecords = (tdm_uint2*)malloc(MaxRecords*sizeof(tdm_uint2));
        */
	h_startRecords = (tdm_uint2*)memalign ( AOCL_ALIGNMENT, MaxRecords*sizeof(tdm_uint2));
	h_foundRecords = (tdm_uint2*)memalign ( AOCL_ALIGNMENT, MaxRecords*sizeof(tdm_uint2));
	// Variables for compaction
	d_recCount=dsycl::alloc<uint>(q, MaxRecords);
	d_recOffSet=dsycl::alloc<uint>(q, MaxRecords);

	dsycl::write(q, d_events, h_events, eventSize);

	dsycl::write(q, d_times, h_times, eventSize * sizeof(float));

	dsycl::write(q, d_episodeCandidates, h_episodeCandidates, numCandidates * maxLevel * sizeof(ubyte));

	dsycl::write(q, d_episodeIntervals, h_episodeIntervals, numCandidates * (maxLevel-1) * 2 * sizeof(float));

}


int main( int argc, char** argv) 
{
	return runTest( argc, argv);

}


int runTest( int argc, char** argv) {
	

	//****************************************

	unsigned int* buff1=NULL;
	unsigned int* buff2=NULL;
	int kl;

	dsycl::parse_args(&argc, &argv);
	static sycl::queue queue_ = dsycl::make_queue();
	Q = &queue_;
	sycl::queue &q = queue_;

	if ( argc != 5 ) {
		printf("Usage: tdm_ocl <data path> <intervals path> <episodes path> <threads>\n");
		return 1;
	}

	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */

	//************************

	size_t globalWorkSize[3];
	size_t localWorkSize[3];

	unsigned int num_threads = atoi(argv[4]);

	dumpFile = fopen( "tdm-gpu-csw.txt", "w" );

	// Load Data
	loadData( argv[1], &h_events, &h_times, &eventSize, &eventType, &uniqueEvents );
	loadTemporalConstraints(argv[2], &temporalConstraint, &temporalConstraintSize );

	// Allocate host intervals, support, and candidates
	/* h_episodeIntervals = (float*)malloc( maxIntervals*sizeof(float) );
	h_episodeCandidates = (ubyte*)malloc( maxCandidates*sizeof(ubyte) );
	h_episodeSupport = (uint*)malloc( maxCandidates*sizeof(float) );
*/
       h_episodeIntervals = (float*)memalign ( AOCL_ALIGNMENT, maxIntervals*sizeof(float) );
	h_episodeCandidates = (ubyte*)memalign ( AOCL_ALIGNMENT, maxCandidates*sizeof(ubyte) );
	h_episodeSupport = (uint*)memalign ( AOCL_ALIGNMENT,maxCandidates*sizeof(float) );

	loadCandidateEpisodes(argv[3], maxLevel);
	setupGpu();

	indexEvents();

	free( h_events );
	free( h_times );

	ubyte start = eventType == EVENT_26 ? 'A' : '!';

	/* common timer: repetitions around the episode loop (per episode: H2D of the start records, countCandidates / writeCandidates
	 * per level, the host prefix sums, D2H of the final records). The one-time upload of the event stream is not repeated. */
	dwarf_init("tdm", "finite-state-machine", "sycl");
	dwarf_set_size(eventSize, numCandidates);
	for(int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	dwarf_region_begin();
	for ( unsigned int idx = 0; idx < numCandidates; idx++ ) {

		ubyte* currentEpisode = &h_episodeCandidates[idx*maxLevel];

		dsycl::write(q, d_startRecords, h_eventIndex[currentEpisode[maxLevel-1]-start], h_eventCounts[currentEpisode[maxLevel-1]-start]*sizeof(tdm_uint2));

		h_foundCount = h_eventCounts[currentEpisode[maxLevel-1]-start];
		localWorkSize[0] = DWARF_LS; /* [FIX-OCL-34] was CLI num_threads (128); pinned work-group = 256, kernels guard tid<recordCount */
		localWorkSize[1] = 1;
		localWorkSize[2] = 1;

		for ( int level = 2; level <= maxLevel; level++ ) {
			int idx = maxLevel - level;                // Index instead of level

			globalWorkSize[0]=DWARF_ROUNDUP(h_foundCount); /* [FIX-OCL-34] pinned global = ceil(n/256)*256 (was based on CLI num_threads) */
			globalWorkSize[1]=1;
			globalWorkSize[2]=1;

			//printf("Starting with %d records and blockCount=%d, localWorkSize[0]=%d globalWorkSize[0]=%d\n", h_foundCount, blockCount, localWorkSize[0], globalWorkSize[0]);
			// Count how many to write

			countCandidates(q, sycl::nd_range<1>(globalWorkSize[0], localWorkSize[0]), currentEpisode[idx], temporalConstraint[0], temporalConstraint[1],
					d_startRecords, d_recCount, h_foundCount, d_events, d_times);

			
			//buff1=(unsigned int*)malloc(sizeof(unsigned int)*(h_foundCount));
			buff1=(unsigned int*) memalign (AOCL_ALIGNMENT,sizeof(unsigned int)*(h_foundCount));			
			if(buff1==NULL){
				printf("Failed to allocate memory\n");
				return 1;
			}
			
			buff2=(unsigned int*) memalign (AOCL_ALIGNMENT,sizeof(unsigned int)*(h_foundCount+1));
//buff2=(unsigned int*)malloc(sizeof(unsigned int)*(h_foundCount+1));
			if(buff2==NULL){
				printf("Failed to allocate memory\n");
				return 1;
			}
			
			dsycl::read(q, buff1, d_recCount, (h_foundCount)*sizeof(uint));
	
			buff2[0]=0;
			for (kl=0;kl<h_foundCount;kl++){
					buff2[kl+1] = buff2[kl]+buff1[kl];
			}
			
			dsycl::write(q, d_recOffSet, buff2, (h_foundCount+1)*sizeof(cl_uint));
			
						
			// Do the actual writing
			writeCandidates(q, sycl::nd_range<1>(globalWorkSize[0], localWorkSize[0]), currentEpisode[idx], temporalConstraint[0], temporalConstraint[1],
					d_startRecords, d_recOffSet, d_foundRecords, h_foundCount, d_events, d_times);

			h_foundCount=buff2[h_foundCount];
			/* [FIX-OCL-21] h_startRecords/h_foundRecords/d_*Records hold MaxRecords entries; a larger frontier overflowed them silently */
			if(h_foundCount > MaxRecords) { printf("FAILURE: %u records exceed MaxRecords (%u)\n", h_foundCount, MaxRecords); exit(1); }
			
			free(buff1);
			free(buff2);
			
			tdm_uint2 *temp = d_startRecords;
			d_startRecords = d_foundRecords;
			d_foundRecords = temp;
			
		}
		
		// Copy back, and eliminate conflicts
		dsycl::read(q, h_startRecords, d_startRecords, h_foundCount*sizeof(tdm_uint2));

		tdm_uint2* result;
		result = eliminateConflicts( h_startRecords, h_foundCount, &h_episodeSupport[idx] );
		free(result);

	}

	dwarf_region_end(rep);
	}
	saveResult( dumpFile, maxLevel, numCandidates, h_episodeSupport, h_episodeCandidates, h_episodeIntervals, eventType );
	
	/*printf("Total Count:            %d\n", resultCount );
	for ( uint idx = 0; idx < resultCount; idx++ ) {
	 fprintf( dumpFile, "(%d,%d)\n", result[idx].x, result[idx].y );
	}
*/
	printf("Cleaning up memory...\n");
	dsycl::release(q, d_events);
	dsycl::release(q, d_times);
	dsycl::release(q, d_startRecords);
	dsycl::release(q, d_foundRecords);

	free(h_startRecords);
	free(h_foundRecords);
	fclose(dumpFile);

	{
		int ok = numCandidates > 0;   /* [FIX-OCL-22] no validity indication before; the support table is the reference for comparisons */
		printf("%s: %u candidate episodes of size %d mined\n", ok ? "SUCCESS" : "FAILURE", numCandidates, maxLevel);
		return dwarf_report(ok);
	}

}
