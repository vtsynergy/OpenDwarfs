/**
 * CUDA translation of the OpenCL port of the CUDA version of
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

#include "dwarf_cuda.h"
#include "dwarf_timer.h"
#include "tdm_kernel.cuh"
#define AOCL_ALIGNMENT 64
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
uint2 *d_startRecords;
uint2  *h_startRecords;

uint2 *d_foundRecords;
uint2  *h_foundRecords;

uint *d_recCount;
uint *d_recOffSet;

uint *d_eventCounts;
uint* h_eventCounts;

uint2** h_eventIndex;

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

	/*h_eventIndex = (uint2**)malloc(uniqueEvents*sizeof(uint2*));
	h_eventCounts = (uint*)malloc(uniqueEvents*sizeof(uint));
        */

        h_eventIndex = (uint2**) memalign ( AOCL_ALIGNMENT, uniqueEvents*sizeof(uint2*));
	h_eventCounts = (uint*) memalign ( AOCL_ALIGNMENT, uniqueEvents*sizeof(uint));

	// First pass to how many events of each type
	memset(h_eventCounts, 0, uniqueEvents*sizeof(unsigned int) );
	for ( long idx = 0; idx < eventSize; idx++ ) {
		h_eventCounts[h_events[idx]-start]++;
	}

	// Generate index arrays
	for ( long idx = 0; idx < uniqueEvents; idx++ ) {
		//printf("%c: %d\n", start+idx, h_eventCounts[idx]);
		//h_eventIndex[idx] = (uint2*)malloc(h_eventCounts[idx]*sizeof(uint2) );
		h_eventIndex[idx] = (uint2*) memalign ( AOCL_ALIGNMENT, h_eventCounts[idx]*sizeof(uint2) );
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
	uint2 i1 = *(uint2*)v1;
	uint2 i2 = *(uint2*)v2;
	return i1.y-i2.y;
}

uint2* eliminateConflicts(uint2* records, uint recordCount, uint* resultCount) {
	//uint2* result = (uint2*)malloc(recordCount*sizeof(uint2));
	uint2* result = (uint2*) memalign ( AOCL_ALIGNMENT,recordCount*sizeof(uint2));
	uint count = 0;
	uint lastRecord = 0;

	// Sort records by endtime
	//qsort( records, recordCount, sizeof(uint2), compareUint2ByY );

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

	CUDA_CHECK(cudaMalloc((void**)&d_events, eventSize * sizeof(ubyte)));
	CUDA_CHECK(cudaMalloc((void**)&d_times, eventSize * sizeof(float)));
	CUDA_CHECK(cudaMalloc((void**)&d_episodeCandidates, maxCandidates* sizeof(ubyte)));
	CUDA_CHECK(cudaMalloc((void**)&d_episodeIntervals, maxIntervals * sizeof(float)));
	CUDA_CHECK(cudaMalloc((void**)&d_episodeSupport, maxCandidates * sizeof(uint)));
	CUDA_CHECK(cudaMalloc((void**)&d_startRecords, MaxRecords*sizeof(uint2)));
	CUDA_CHECK(cudaMalloc((void**)&d_foundRecords, MaxRecords*sizeof(uint2)));
	/*
	h_startRecords = (uint2*)malloc(MaxRecords*sizeof(uint2));
	h_foundRecords = (uint2*)malloc(MaxRecords*sizeof(uint2));
        */
	h_startRecords = (uint2*)memalign ( AOCL_ALIGNMENT, MaxRecords*sizeof(uint2));
	h_foundRecords = (uint2*)memalign ( AOCL_ALIGNMENT, MaxRecords*sizeof(uint2));
	// Variables for compaction
	CUDA_CHECK(cudaMalloc((void**)&d_recCount, MaxRecords*sizeof(uint)));
	CUDA_CHECK(cudaMalloc((void**)&d_recOffSet, MaxRecords*sizeof(uint)));

#define TDM_COPY(dst, src, bytes, kind) do { CUDA_CHECK(cudaMemcpy(dst, src, bytes, kind)); CUDA_CHECK(cudaDeviceSynchronize()); } while (0)
	TDM_COPY(d_events, h_events, eventSize, cudaMemcpyHostToDevice);
	TDM_COPY(d_times, h_times, eventSize * sizeof(float), cudaMemcpyHostToDevice);
	TDM_COPY(d_episodeCandidates, h_episodeCandidates, numCandidates * maxLevel * sizeof(ubyte), cudaMemcpyHostToDevice);
	TDM_COPY(d_episodeIntervals, h_episodeIntervals, numCandidates * (maxLevel-1) * 2 * sizeof(float), cudaMemcpyHostToDevice);

}


int main( int argc, char** argv) 
{
	return runTest( argc, argv);

}


int runTest( int argc, char** argv) {
	

	//****************************************

	int kl;
	
	unsigned int* buff1=NULL;
	unsigned int* buff2=NULL;

	dwarf_cuda_init(&argc, &argv);
	
	if ( argc != 5 ) {
		printf("Usage: tdm_ocl <data path> <intervals path> <episodes path> <threads>\n");
		return 1;
	}

	//************************


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
	dwarf_init("tdm", "finite-state-machine", "cuda");
	dwarf_set_size(eventSize, numCandidates);
	for(int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	dwarf_region_begin();
	for ( unsigned int idx = 0; idx < numCandidates; idx++ ) {

		ubyte* currentEpisode = &h_episodeCandidates[idx*maxLevel];

		TDM_COPY(d_startRecords, h_eventIndex[currentEpisode[maxLevel-1]-start], h_eventCounts[currentEpisode[maxLevel-1]-start]*sizeof(uint2), cudaMemcpyHostToDevice);

		h_foundCount = h_eventCounts[currentEpisode[maxLevel-1]-start];
		/* free-size kernels: DWARF_LS threads per block (the CLI thread count is ignored, as in OpenCL); guard tid<recordCount */

		for ( int level = 2; level <= maxLevel; level++ ) {
			int idx = maxLevel - level;                // Index instead of level

			// Count how many to write
			countCandidates<<<DWARF_NBLOCKS(h_foundCount), DWARF_LS>>>(currentEpisode[idx], temporalConstraint[0], temporalConstraint[1],
					d_startRecords, d_recCount, h_foundCount, d_events, d_times);
			CUDA_SYNC_LAUNCH();

			
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
			
			TDM_COPY(buff1, d_recCount, (h_foundCount)*sizeof(uint), cudaMemcpyDeviceToHost);
	
			buff2[0]=0;
			for (kl=0;kl<h_foundCount;kl++){
					buff2[kl+1] = buff2[kl]+buff1[kl];
			}
			
			TDM_COPY(d_recOffSet, buff2, (h_foundCount+1)*sizeof(uint), cudaMemcpyHostToDevice);
			
						
			// Do the actual writing
			writeCandidates<<<DWARF_NBLOCKS(h_foundCount), DWARF_LS>>>(currentEpisode[idx], temporalConstraint[0], temporalConstraint[1],
					d_startRecords, d_recOffSet, d_foundRecords, h_foundCount, d_events, d_times);
			CUDA_SYNC_LAUNCH();

			h_foundCount=buff2[h_foundCount];
			/* [FIX-OCL-21] h_startRecords/h_foundRecords/d_*Records hold MaxRecords entries; a larger frontier overflowed them silently */
			if(h_foundCount > MaxRecords) { printf("FAILURE: %u records exceed MaxRecords (%u)\n", h_foundCount, MaxRecords); exit(1); }
			
			free(buff1);
			free(buff2);
			
			uint2 *temp = d_startRecords;
			d_startRecords = d_foundRecords;
			d_foundRecords = temp;
			
		}
		
		// Copy back, and eliminate conflicts
		TDM_COPY(h_startRecords, d_startRecords, h_foundCount*sizeof(uint2), cudaMemcpyDeviceToHost);

		uint2* result;
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
	cudaFree(d_events);
	cudaFree(d_times);
	cudaFree(d_startRecords);
	cudaFree(d_foundRecords);

	free(h_startRecords);
	free(h_foundRecords);
	fclose(dumpFile);


	{
		int ok = numCandidates > 0;   /* [FIX-OCL-22] no validity indication before; the support table is the reference for comparisons */
		printf("%s: %u candidate episodes of size %d mined\n", ok ? "SUCCESS" : "FAILURE", numCandidates, maxLevel);
		return dwarf_report(ok);
	}

}
