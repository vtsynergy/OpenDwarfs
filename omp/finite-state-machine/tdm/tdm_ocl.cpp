/**
 * Temporal Data Mining (finite-state-machine dwarf) -- OpenMP target-offload translation of
 * opencl/finite-state-machine/tdm/tdm_ocl.cpp + tdm_ocl_kernel.cl (ported to OpenCL from the CUDA GPU
 * Temporal Data Mining by Sean Ponce).
 *
 * Same flow: read the event stream, temporal constraints and candidate episodes; index the events by type; for each
 * candidate episode walk the levels backwards: countCandidates (kernel) -> prefix sum of the counts on the host ->
 * writeCandidates (kernel); the found records become the next level's start records. The support of an episode is
 * the number of non-conflicting (non-overlapping) final records (eliminateConflicts, host).
 * Usage: tdm <events.csv> <intervals.txt> <episodes.txt> <threads>   -> writes tdm-gpu-csw.txt
 * Timed region: everything from the device uploads to the last episode (the device loop with its small host steps).
 * Validation (added): the support counts in tdm-gpu-csw.txt are diffed against the OpenCL run; the app itself
 * checks that every episode produced a support.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>
#include "types.h"
#include "global.h"
#include "dataio.h"
#include "dwarf_timer.h"
#include "dwarf_launch.h" /* [FIX-OMP-13] pinned launch config */

typedef struct { unsigned int x, y; } uint2;

unsigned int numCandidates;
int episodesCulled = 0;
ubyte* h_events;
float* h_times;
int maxLevel;
ubyte* h_episodeCandidates;
float* h_episodeIntervals;
uint* h_episodeSupport;
uint eventSize;
uint eventType, uniqueEvents;
unsigned int h_foundCount;
uint2 *h_startRecords, *h_foundRecords;
uint* h_eventCounts;
uint2** h_eventIndex;
FILE* dumpFile;
float* temporalConstraint;
unsigned int temporalConstraintSize;

/* "device" buffers: host arrays that carry the OpenMP mapping (d_* of the OpenCL code) */
static ubyte* d_events;
static float* d_times;
static uint2 *d_startRecords, *d_foundRecords;
static uint *d_recCount, *d_recOffSet;

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

        h_eventIndex = (uint2**) malloc(uniqueEvents*sizeof(uint2*));
	h_eventCounts = (uint*) malloc(uniqueEvents*sizeof(uint));

	// First pass to how many events of each type
	memset(h_eventCounts, 0, uniqueEvents*sizeof(unsigned int) );
	for ( long idx = 0; idx < eventSize; idx++ ) {
		h_eventCounts[h_events[idx]-start]++;
	}

	// Generate index arrays
	for ( long idx = 0; idx < uniqueEvents; idx++ ) {
		//printf("%c: %d\n", start+idx, h_eventCounts[idx]);
		//h_eventIndex[idx] = (uint2*)malloc(h_eventCounts[idx]*sizeof(uint2) );
		h_eventIndex[idx] = (uint2*) malloc(h_eventCounts[idx]*sizeof(uint2) );
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
	uint2* result = (uint2*) malloc(recordCount*sizeof(uint2));
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


/* ---- kernel countCandidates ---- */
static void k_count(ubyte targetEvent, float minInterval, float maxInterval, uint2* startRecords, uint* foundRecordsCount,
                    uint recordCount, ubyte* ev, float* tm, size_t nev) {
/* [FIX-OMP-13] free-size: pinned 256 x ceil(recordCount/256) = OpenCL */
#pragma omp target teams distribute parallel for DWARF_LAUNCH(recordCount) map(present, alloc : startRecords[0 : MaxRecords], foundRecordsCount[0 : MaxRecords], ev[0 : nev], tm[0 : nev])
  for (uint tid = 0; tid < recordCount; tid++) {
    uint2 ep;  // .x is start, .y is end
    uint myIdx;
    float startTime, curTime;
    unsigned int listSize = 0;
    // Retrieve endIndex from startRecords and set to startIdx
    ep = startRecords[tid];
    myIdx = ep.x;
    // Retrieve interval for this event
    curTime = startTime = tm[myIdx];
    // Advance through events until minimum interval has been reached
    while (startTime - curTime < minInterval && myIdx > 0) {
      myIdx--;
      curTime = tm[myIdx];
    }
    // While the current event has not exceeded the maximum interval
    while (startTime - curTime <= maxInterval && myIdx > 0) {
      if (ev[myIdx] == targetEvent) listSize++;
      myIdx--;
      curTime = tm[myIdx];
    }
    foundRecordsCount[tid] = listSize;
  }
}

/* ---- kernel writeCandidates ---- */
static void k_write(ubyte targetEvent, float minInterval, float maxInterval, uint2* startRecords, uint* startOffset,
                    uint2* foundRecords, uint recordCount, ubyte* ev, float* tm, size_t nev) {
/* [FIX-OMP-13] free-size: pinned 256 x ceil(recordCount/256) = OpenCL */
#pragma omp target teams distribute parallel for DWARF_LAUNCH(recordCount) map(present, alloc : startRecords[0 : MaxRecords], startOffset[0 : MaxRecords], foundRecords[0 : MaxRecords], ev[0 : nev], tm[0 : nev])
  for (uint tid = 0; tid < recordCount; tid++) {
    uint2 ep;
    uint myIdx, offset;
    float startTime, curTime;
    unsigned int listSize = 0;
    ep = startRecords[tid];
    offset = startOffset[tid];
    myIdx = ep.x;
    curTime = startTime = tm[myIdx];
    while (startTime - curTime < minInterval && myIdx > 0) {
      myIdx--;
      curTime = tm[myIdx];
    }
    while (startTime - curTime <= maxInterval && myIdx > 0) {
      if (ev[myIdx] == targetEvent) {
        foundRecords[offset + listSize].x = myIdx;
        foundRecords[offset + listSize].y = ep.y;
        listSize++;
      }
      myIdx--;
      curTime = tm[myIdx];
    }
  }
}

int main(int argc, char** argv) {
  unsigned int* buff1 = NULL;
  unsigned int* buff2 = NULL;
  int kl;
  if (argc != 5) {
    printf("Usage: tdm_ocl <data path> <intervals path> <episodes path> <threads>\n");
    return 1;
  }
  unsigned int num_threads = atoi(argv[4]);
  (void)num_threads;  // work-group size of the OpenCL launch; irrelevant for the flat parallel loops
  dumpFile = fopen("tdm-gpu-csw.txt", "w");
  loadData(argv[1], &h_events, &h_times, &eventSize, &eventType, &uniqueEvents);
  loadTemporalConstraints(argv[2], &temporalConstraint, &temporalConstraintSize);
  h_episodeIntervals = (float*)malloc(maxIntervals * sizeof(float));
  h_episodeCandidates = (ubyte*)malloc(maxCandidates * sizeof(ubyte));
  h_episodeSupport = (uint*)malloc(maxCandidates * sizeof(float));
  loadCandidateEpisodes(argv[3], maxLevel);

  /* setupGpu(): host mirrors of the device buffers */
  d_events = (ubyte*)malloc(eventSize);
  d_times = (float*)malloc(eventSize * sizeof(float));
  d_startRecords = (uint2*)malloc(MaxRecords * sizeof(uint2));
  d_foundRecords = (uint2*)malloc(MaxRecords * sizeof(uint2));
  d_recCount = (uint*)malloc(MaxRecords * sizeof(uint));
  d_recOffSet = (uint*)malloc(MaxRecords * sizeof(uint));
  h_startRecords = (uint2*)malloc(MaxRecords * sizeof(uint2));
  h_foundRecords = (uint2*)malloc(MaxRecords * sizeof(uint2));
  memcpy(d_events, h_events, eventSize);
  memcpy(d_times, h_times, eventSize * sizeof(float));
  indexEvents();
  free(h_events);
  free(h_times);

  ubyte start = eventType == EVENT_26 ? 'A' : '!';
  size_t nev = eventSize;
  ubyte* dev = d_events; float* dtm = d_times; uint2 *dstart = d_startRecords, *dfound = d_foundRecords;
  uint *dcount = d_recCount, *doff = d_recOffSet;

  dwarf_init("tdm", "finite-state-machine", "omp");
  dwarf_set_size(eventSize, numCandidates);
  /* [FIX-OMP-14] OpenCL creates its buffers before the timed region and only transfers inside it: allocate once outside the region, `target update` inside */
  /* events and times are uploaded once, as in the OpenCL setupGpu() */
#pragma omp target enter data map(to : dev[0 : nev], dtm[0 : nev]) map(alloc : dstart[0 : MaxRecords], dfound[0 : MaxRecords], dcount[0 : MaxRecords], doff[0 : MaxRecords])
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    d_startRecords = dstart; d_foundRecords = dfound;  // undo the pointer swaps of the previous repetition
    dwarf_region_begin();
    for (unsigned int idx = 0; idx < numCandidates; idx++) {
      ubyte* currentEpisode = &h_episodeCandidates[idx * maxLevel];
      unsigned int n0 = h_eventCounts[currentEpisode[maxLevel - 1] - start];
      uint2* sr = d_startRecords;
      memcpy(sr, h_eventIndex[currentEpisode[maxLevel - 1] - start], n0 * sizeof(uint2));
#pragma omp target update to(sr[0 : n0])
      h_foundCount = n0;
      for (int level = 2; level <= maxLevel; level++) {
        int eidx = maxLevel - level;  // Index instead of level
        uint2* srec = d_startRecords; uint2* frec = d_foundRecords;
        k_count(currentEpisode[eidx], temporalConstraint[0], temporalConstraint[1], srec, dcount, h_foundCount, dev, dtm, nev);
        buff1 = (unsigned int*)malloc(sizeof(unsigned int) * (h_foundCount + 1));
        buff2 = (unsigned int*)malloc(sizeof(unsigned int) * (h_foundCount + 1));
        if (buff1 == NULL || buff2 == NULL) { printf("Failed to allocate memory\n"); return 1; }
        unsigned int fc = h_foundCount;
        // copy counts to the host, prefix sum, copy offsets back (as the OpenCL host code)
#pragma omp target update from(dcount[0 : fc])
        memcpy(buff1, dcount, fc * sizeof(uint));
        buff2[0] = 0;
        for (kl = 0; kl < (int)fc; kl++) buff2[kl + 1] = buff2[kl] + buff1[kl];
        memcpy(doff, buff2, (fc + 1) * sizeof(uint));
#pragma omp target update to(doff[0 : fc + 1])
        k_write(currentEpisode[eidx], temporalConstraint[0], temporalConstraint[1], srec, doff, frec, h_foundCount, dev, dtm, nev);
        h_foundCount = buff2[h_foundCount];
        /* [FIX-OMP-03] the record buffers hold MaxRecords entries; a larger frontier overflowed them silently (FIX-OCL-21) */
        if (h_foundCount > MaxRecords) { printf("FAILURE: %u records exceed MaxRecords (%u)\n", h_foundCount, MaxRecords); return 1; }
        free(buff1);
        free(buff2);
        uint2* temp = d_startRecords;
        d_startRecords = d_foundRecords;
        d_foundRecords = temp;
      }
      {
        uint2* sr2 = d_startRecords; unsigned int fc2 = h_foundCount;
#pragma omp target update from(sr2[0 : fc2])
        memcpy(h_startRecords, sr2, fc2 * sizeof(uint2));
      }
      uint2* result = eliminateConflicts(h_startRecords, h_foundCount, &h_episodeSupport[idx]);
      free(result);
    }
    dwarf_region_end(rep);
  }
#pragma omp target exit data map(delete : dev[0 : nev], dtm[0 : nev], dstart[0 : MaxRecords], dfound[0 : MaxRecords], dcount[0 : MaxRecords], doff[0 : MaxRecords])
  saveResult(dumpFile, maxLevel, numCandidates, h_episodeSupport, h_episodeCandidates, h_episodeIntervals, eventType);
  printf("Cleaning up memory...\n");
  fclose(dumpFile);
  int ok = numCandidates > 0;
  printf("%s: %u candidate episodes of size %d mined\n", ok ? "SUCCESS" : "FAILURE", numCandidates, maxLevel);
  return dwarf_report(ok);
}
