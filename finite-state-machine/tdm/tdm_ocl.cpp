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

#include "types.h"
#include "global.h"
#include "dataio.h"

#include "../../include/rdtsc.h"
#include "../../include/common_args.h"

unsigned int numCandidates;
int episodesCulled = 0;

ubyte* h_events;
cl_mem d_events;

float* h_times;
cl_mem d_times;

int maxLevel;

ubyte* h_episodeCandidates;
cl_mem d_episodeCandidates;

float* h_episodeIntervals;
cl_mem d_episodeIntervals;

uint* h_episodeSupport;
cl_mem d_episodeSupport;

int* h_episodeIndices;
cl_mem d_episodeIndices;

uint eventSize;
uint eventType, uniqueEvents;

unsigned int  h_foundCount;
cl_mem d_startRecords;
cl_uint2  *h_startRecords;

cl_mem d_foundRecords;
cl_uint2  *h_foundRecords;

cl_mem d_recCount;
cl_mem d_recOffSet;

cl_mem d_eventCounts;
uint* h_eventCounts;

cl_mem d_eventIndex;
cl_uint2** h_eventIndex;

FILE* dumpFile;

float* temporalConstraint;
unsigned int temporalConstraintSize;

void runTest( int argc, char** argv);

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

	h_eventIndex = (cl_uint2**)malloc(uniqueEvents*sizeof(cl_uint2*));
	h_eventCounts = (uint*)malloc(uniqueEvents*sizeof(uint));

	// First pass to how many events of each type
	memset(h_eventCounts, 0, uniqueEvents*sizeof(unsigned int) );
	for ( long idx = 0; idx < eventSize; idx++ ) {
		h_eventCounts[h_events[idx]-start]++;
	}

	// Generate index arrays
	for ( long idx = 0; idx < uniqueEvents; idx++ ) {
		//printf("%c: %d\n", start+idx, h_eventCounts[idx]);
		h_eventIndex[idx] = (cl_uint2*)malloc(h_eventCounts[idx]*sizeof(cl_uint2) );
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
	cl_uint2 i1 = *(cl_uint2*)v1;
	cl_uint2 i2 = *(cl_uint2*)v2;
	return i1.y-i2.y;
}

cl_uint2* eliminateConflicts(cl_uint2* records, uint recordCount, uint* resultCount) {
	cl_uint2* result = (cl_uint2*)malloc(recordCount*sizeof(cl_uint2));
	uint count = 0;
	uint lastRecord = 0;

	// Sort records by endtime
	//qsort( records, recordCount, sizeof(cl_uint2), compareUint2ByY );

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

	cl_int errcode;

	d_events=clCreateBuffer(context, CL_MEM_READ_WRITE, eventSize * sizeof(ubyte), NULL, &errcode);
	d_times=clCreateBuffer(context, CL_MEM_READ_WRITE, eventSize * sizeof(float), NULL, &errcode);
	d_episodeCandidates=clCreateBuffer(context, CL_MEM_READ_WRITE, maxCandidates* sizeof(ubyte), NULL, &errcode);
	d_episodeIntervals=clCreateBuffer(context, CL_MEM_READ_WRITE, maxIntervals * sizeof(float), NULL, &errcode);
	d_episodeSupport=clCreateBuffer(context, CL_MEM_READ_WRITE, maxCandidates * sizeof(uint), NULL, &errcode);
	d_startRecords=clCreateBuffer(context, CL_MEM_READ_WRITE, MaxRecords*sizeof(cl_uint2), NULL, &errcode);
	d_foundRecords=clCreateBuffer(context, CL_MEM_READ_WRITE, MaxRecords*sizeof(cl_uint2), NULL, &errcode);
	h_startRecords = (cl_uint2*)malloc(MaxRecords*sizeof(cl_uint2));
	h_foundRecords = (cl_uint2*)malloc(MaxRecords*sizeof(cl_uint2));

	// Variables for compaction
	d_recCount=clCreateBuffer(context, CL_MEM_READ_WRITE, MaxRecords*sizeof(cl_uint), NULL, &errcode);
	d_recOffSet=clCreateBuffer(context, CL_MEM_READ_WRITE, MaxRecords*sizeof(cl_uint), NULL, &errcode);

	errcode = clEnqueueWriteBuffer(commands, d_events, CL_TRUE, 0, eventSize, (void *) h_events, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "d_events Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue write buffer!");

	errcode = clEnqueueWriteBuffer(commands, d_times, CL_TRUE, 0, eventSize * sizeof(float), (void *) h_times, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "d_times Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue write buffer!");	

	errcode = clEnqueueWriteBuffer(commands, d_episodeCandidates, CL_TRUE, 0, numCandidates * maxLevel * sizeof(ubyte), (void *) h_episodeCandidates, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "d_episodeCandidates Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue write buffer!");

	errcode = clEnqueueWriteBuffer(commands, d_episodeIntervals, CL_TRUE, 0, numCandidates * (maxLevel-1) * 2 * sizeof(float), (void *) h_episodeIntervals, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "d_episodeIntervals Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue write buffer!");

}


int main( int argc, char** argv) 
{
	runTest( argc, argv);

}


void runTest( int argc, char** argv) {
	

	//****************************************

	cl_program clProgram;
	cl_kernel clKernel_writeCandidates;
	cl_kernel clKernel_countCandidates;
	cl_int dev_type;
	cl_int errcode;
	FILE *kernelFile;
	char *kernelSource;
	size_t kernelLength;
	int kl;
	
	unsigned int* buff1=NULL;
	unsigned int* buff2=NULL;

	ocd_init(&argc, &argv, NULL);
	ocd_initCL();
	
	if ( argc != 5 ) {
		printf("Usage: tdm_ocl <data path> <intervals path> <episodes path> <threads>\n");
		return;
	}

	kernelFile = fopen("tdm_ocl_kernel.cl", "r");
	fseek(kernelFile, 0, SEEK_END);
	kernelLength = (size_t) ftell(kernelFile);
	kernelSource = (char *) malloc(sizeof(char)*kernelLength);
	rewind(kernelFile);
	fread((void *) kernelSource, kernelLength, 1, kernelFile);
	fclose(kernelFile);

	clProgram = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, &kernelLength, &errcode);
	CHKERR(errcode, "Failed to create program with source!");

	free(kernelSource);	
	errcode = clBuildProgram(clProgram, 1, &device_id, NULL, NULL, NULL);
	if (errcode == CL_BUILD_PROGRAM_FAILURE)                                                                                                                                       
	{
		char *log;
		size_t logLength;
		errcode = clGetProgramBuildInfo(clProgram, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLength);
		log = (char *) malloc(sizeof(char)*logLength);
		errcode = clGetProgramBuildInfo(clProgram, device_id, CL_PROGRAM_BUILD_LOG, logLength, (void *) log, NULL);
		fprintf(stderr, "Kernel build error! Log:\n%s", log);
		free(log);
		return;
	}
	CHKERR(errcode, "Failed to build program!");

	clKernel_writeCandidates = clCreateKernel(clProgram, "writeCandidates", &errcode);
	CHKERR(errcode, "Failed to create kernel!");
	clKernel_countCandidates = clCreateKernel(clProgram, "countCandidates", &errcode);
	CHKERR(errcode, "Failed to create kernel!");

	//************************

	size_t globalWorkSize[3];
	size_t localWorkSize[3];

	unsigned int num_threads = atoi(argv[4]);

	dumpFile = fopen( "tdm-gpu-csw.txt", "w" );

	// Load Data
	loadData( argv[1], &h_events, &h_times, &eventSize, &eventType, &uniqueEvents );
	loadTemporalConstraints(argv[2], &temporalConstraint, &temporalConstraintSize );

	// Allocate host intervals, support, and candidates
	h_episodeIntervals = (float*)malloc( maxIntervals*sizeof(float) );
	h_episodeCandidates = (ubyte*)malloc( maxCandidates*sizeof(ubyte) );
	h_episodeSupport = (uint*)malloc( maxCandidates*sizeof(float) );

	loadCandidateEpisodes(argv[3], maxLevel);
	setupGpu();

	indexEvents();

	free( h_events );
	free( h_times );

	ubyte start = eventType == EVENT_26 ? 'A' : '!';

	for ( unsigned int idx = 0; idx < numCandidates; idx++ ) {

		ubyte* currentEpisode = &h_episodeCandidates[idx*maxLevel];

		errcode = clEnqueueWriteBuffer(commands, d_startRecords, CL_TRUE, 0, h_eventCounts[currentEpisode[maxLevel-1]-start]*sizeof(cl_uint2), (void *) h_eventIndex[currentEpisode[maxLevel-1]-start], 0, NULL, &ocdTempEvent);
		clFinish(commands);
		START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "d_startRecords Copy", ocdTempTimer)
		END_TIMER(ocdTempTimer)
		CHKERR(errcode, "Failed to enqueue write buffer!");

		h_foundCount = h_eventCounts[currentEpisode[maxLevel-1]-start];
		localWorkSize[0] = num_threads;
		localWorkSize[1] = 1;
		localWorkSize[2] = 1;

		for ( int level = 2; level <= maxLevel; level++ ) {
			int idx = maxLevel - level;                // Index instead of level

			uint blockCount = (h_foundCount/num_threads)+(h_foundCount % num_threads == 0 ? 0 : 1);
			globalWorkSize[0]=blockCount*num_threads;
			globalWorkSize[1]=1;
			globalWorkSize[2]=1;

			//printf("Starting with %d records and blockCount=%d, localWorkSize[0]=%d globalWorkSize[0]=%d\n", h_foundCount, blockCount, localWorkSize[0], globalWorkSize[0]);
			// Count how many to write

			errcode = clSetKernelArg(clKernel_countCandidates, 0, sizeof(ubyte), (void *) &currentEpisode[idx]);
			errcode |= clSetKernelArg(clKernel_countCandidates, 1, sizeof(float), (void *) &temporalConstraint[0]);
			errcode |= clSetKernelArg(clKernel_countCandidates, 2, sizeof(float), (void *) &temporalConstraint[1]);
			errcode = clSetKernelArg(clKernel_countCandidates, 3, sizeof(cl_uint2*), (void *) &d_startRecords);
			errcode |= clSetKernelArg(clKernel_countCandidates, 4, sizeof(uint*), (void *) &d_recCount);
			errcode |= clSetKernelArg(clKernel_countCandidates, 5, sizeof(uint), (void *) &h_foundCount);
			errcode |= clSetKernelArg(clKernel_countCandidates, 6, sizeof(ubyte*), (void *) &d_events);
			errcode |= clSetKernelArg(clKernel_countCandidates, 7, sizeof(float*), (void *) &d_times);
			CHKERR(errcode, "Failed to set kernel arguments (1)!");

			errcode = clEnqueueNDRangeKernel(commands, clKernel_countCandidates, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
			clFinish(commands);
			CHKERR(errcode, "Failed to enqueue kernel clKernel_countCandidates!");			
			START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "countCandidates kernel", ocdTempTimer)
			END_TIMER(ocdTempTimer)

			
			buff1=(unsigned int*)malloc(sizeof(unsigned int)*(h_foundCount));
			if(buff1==NULL){
				printf("Failed to allocate memory\n");
				return;
			}
			
			buff2=(unsigned int*)malloc(sizeof(unsigned int)*(h_foundCount+1));
			if(buff2==NULL){
				printf("Failed to allocate memory\n");
				return;
			}
			
			errcode = clEnqueueReadBuffer(commands, d_recCount, CL_TRUE, 0, (h_foundCount)*sizeof(uint), (void *) buff1, 0, NULL, &ocdTempEvent);
			clFinish(commands);
			START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "d_recCount Copy", ocdTempTimer)
			END_TIMER(ocdTempTimer)
			CHKERR(errcode, "Failed to enqueue read buffer!");
	
			buff2[0]=0;
			for (kl=0;kl<h_foundCount;kl++){
					buff2[kl+1] = buff2[kl]+buff1[kl];
			}
			
			errcode = clEnqueueWriteBuffer(commands, d_recOffSet, CL_TRUE, 0, (h_foundCount+1)*sizeof(cl_uint), (void *) buff2, 0, NULL, &ocdTempEvent);
			clFinish(commands);
			START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "d_recOffset Copy", ocdTempTimer)
			END_TIMER(ocdTempTimer)
			CHKERR(errcode, "Failed to enqueue write buffer!");
			
						
			// Do the actual writing
			errcode = clSetKernelArg(clKernel_writeCandidates, 0, sizeof(ubyte), (void *) &currentEpisode[idx]);
			errcode |= clSetKernelArg(clKernel_writeCandidates, 1, sizeof(float), (void *) &temporalConstraint[0]);
			errcode |= clSetKernelArg(clKernel_writeCandidates, 2, sizeof(float), (void *) &temporalConstraint[1]);
			errcode = clSetKernelArg(clKernel_writeCandidates, 3, sizeof(cl_uint2*), (void *) &d_startRecords);
			errcode |= clSetKernelArg(clKernel_writeCandidates, 4, sizeof(unsigned int*), (void *) &d_recOffSet);
			errcode |= clSetKernelArg(clKernel_writeCandidates, 5, sizeof(cl_uint2*), (void *) &d_foundRecords);
			errcode |= clSetKernelArg(clKernel_writeCandidates, 6, sizeof(uint), (void *) &h_foundCount);
			errcode |= clSetKernelArg(clKernel_writeCandidates, 7, sizeof(ubyte*), (void *) &d_events);
			errcode |= clSetKernelArg(clKernel_writeCandidates, 8, sizeof(float*), (void *) &d_times);
			CHKERR(errcode, "Failed to set kernel arguments (2)!");

			errcode = clEnqueueNDRangeKernel(commands, clKernel_writeCandidates, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
			clFinish(commands);
			CHKERR(errcode, "Failed to enqueue kernel clKernel_writeCandidates!");
			START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "writeCandidates kernel", ocdTempTimer)
			END_TIMER(ocdTempTimer)

			h_foundCount=buff2[h_foundCount];
			
			free(buff1);
			free(buff2);
			
			cl_mem temp = d_startRecords;
			d_startRecords = d_foundRecords;
			d_foundRecords = temp;
			
		}
		
		// Copy back, and eliminate conflicts
		errcode = clEnqueueReadBuffer(commands, d_startRecords, CL_TRUE, 0, h_foundCount*sizeof(cl_uint2), (void *) h_startRecords, 0, NULL, &ocdTempEvent);
		clFinish(commands);
		START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "d_startRecords copy", ocdTempTimer)
		END_TIMER(ocdTempTimer)

		cl_uint2* result;
		result = eliminateConflicts( h_startRecords, h_foundCount, &h_episodeSupport[idx] );
		free(result);

	}

	saveResult( dumpFile, maxLevel, numCandidates, h_episodeSupport, h_episodeCandidates, h_episodeIntervals, eventType );
	
	//printf("Total Count:            %d\n", resultCount );
	//for ( uint idx = 0; idx < resultCount; idx++ ) {
	//  fprintf( dumpFile, "(%d,%d)\n", result[idx].x, result[idx].y );
	//}

	//printf("Cleaning up memory...\n");
	clReleaseMemObject(d_events);
	clReleaseMemObject(d_times);
	clReleaseMemObject(d_startRecords);
	clReleaseMemObject(d_foundRecords);

	free(h_startRecords);
	free(h_foundRecords);
	fclose(dumpFile);

	clReleaseKernel(clKernel_countCandidates);
	clReleaseKernel(clKernel_writeCandidates);
	clReleaseProgram(clProgram);
	clReleaseCommandQueue(commands);
	clReleaseContext(context);

	ocd_finalize();

}
