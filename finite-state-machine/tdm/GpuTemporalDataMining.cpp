/*
 * GPU Temporal Data Mining
 * Author: Sean Ponce
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "../../include/rdtsc.h"

int platform_id = PLATFORM_ID, n_device = DEVICE_ID;
// includes, project


#define CHKERR(err, str) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "CL Error %d: %s\n", err, str); \
        exit(1); \
    }

//#define USEGPU 1

#define min(x,y) (x < y ? x : y)

// includes, kernels
//#include <GpuTemporalDataMining_kernel.h>

#include "types.h"
#include "global.h"
#include "OptimalThresholds.h"


unsigned int numCandidates;
int episodesCulled = 0;

UBYTE* h_events;
cl_mem d_events;

float* h_times;
cl_mem d_times;

UBYTE* h_episodeCandidates;
cl_mem d_episodeCandidates;

UBYTE* h_episodeCandidatesBuffer;
cl_mem d_episodeCandidatesBuffer;

float* h_episodeIntervals;
cl_mem d_episodeIntervals;

float* h_episodeIntervalsBuffer;
cl_mem d_episodeIntervalsBuffer;

unsigned int candidateWidth, candidateHeight;
unsigned int eventWidth, eventHeight;

float* h_mapRecords;	// Debugging only
cl_mem d_mapRecords;

float* h_episodeSupport;
cl_mem d_episodeSupport;
bool* h_episodePairs;
cl_mem d_episodePairs;
int* h_episodeIndices;
cl_mem d_episodeIndices;

long eventSize;
long padEventSize;

unsigned int timer, generating_timer, a1_counting_timer, a2_counting_timer, total_timer;

FILE* dumpFile;

float* temporalConstraint;
unsigned int temporalConstraintSize;

int symbolSize;

size_t MaxImageWidth;
size_t MaxImageHeight;

// OpenCL variables
cl_device_id     device_id;
cl_context       context;
cl_command_queue commands;
cl_program       program;
cl_kernel        kernel_countCandidates;
cl_kernel        kernel_countCandidatesStatic;
cl_kernel        kernel_countCandidatesMapMerge;
cl_kernel        kernel_countCandidatesMapMergeStatic;
cl_kernel        kernel_analyzeSupport;
cl_kernel        kernel_generateEpisodeCandidatesKernel;

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest( int argc, char** argv);
void getDeviceVariables(cl_device_id device);

extern "C"
void computeGold( float* reference, float* idata, const unsigned int len);

char symbolToChar( char l, int n )
{
	return (l - 'a')*8 + (n-1) + '!';
}

void charToSymbol( char c, char& l, int& n )
{
	n = (c-'!')%8+1;
	l = (c-'!'-(n-1))/8 + 'a';

}

size_t roundup(size_t size)
{
	size_t width = size <= MaxImageWidth ? size : MaxImageWidth;
	size_t height = (size + MaxImageWidth - 1) / MaxImageWidth;
	return width * height;
}

void initEpisodeCandidates()
{
	if ( symbolSize == ONE_CHAR )
	{
		numCandidates = uniqueEvents;
		for ( char i = 'A'; i <= 'Z'; i++ )
		{
			h_episodeCandidates[i-'A'] = i;
		}
	}
	else
	{
		numCandidates = 64;
		for ( unsigned char i = 0; i < 64; i++ )
		{
			h_episodeCandidates[i] = i + '!';
		}
	}
}


int bindTexture(int offset, cl_mem* texture, cl_mem memory, size_t size, cl_channel_type dataType)
{
    size_t origin[3];
    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 0;
    size_t region[3];
    region[0] = size <= MaxImageWidth ? size : MaxImageWidth;
    region[1] = (size + MaxImageWidth - 1) / MaxImageWidth;
    region[2] = 1;

    if(region[0] == 0)
        region[0] = 1;
    if(region[1] == 0)
        region[1] = 1;

    cl_image_format format;
    format.image_channel_order = CL_R;
    format.image_channel_data_type = dataType;

    int err = CL_SUCCESS;
    *texture = clCreateImage2D(context, CL_MEM_READ_ONLY, &format, region[0], region[1], 0, NULL, &err);//region[0], region[1], 0, NULL, &err);
    CHKERR(err, "Unable to create texture!");

    if(size != 0)
	    err = clEnqueueCopyBufferToImage(commands, memory, *texture, 0, origin, region, 0, NULL, NULL);
    CHKERR(err, "Unable to buffer texture!");
    return CL_SUCCESS;
}

int unbindTexture(cl_mem* texture, cl_mem memory, size_t size)
{
//    size_t origin[3];
//    origin[0] = 0;
//    origin[1] = 0;
//    origin[2] = 0;
//    size_t region[3];
//    region[0] = size <= MaxImageWidth ? size : MaxImageWidth;
//    region[1] = (size + MaxImageWidth - 1) / MaxImageWidth;
//    region[2] = 0;
//
//    if(region[0] == 0)
//        region[0] = 1;
//    if(region[1] == 0)
//        region[1] = 1;
//    int err = clEnqueueCopyImageToBuffer(commands, *texture, memory, origin, region, 0, 0, NULL, NULL);
//    CHKERR(err, "Unable to copy image to buffer!");
    return clReleaseMemObject(*texture);
}

int getNextValidCandidate(int prefixLength, int currentIdx, int nextIdx)
{
	//nextIdx++;
	for ( int idx = nextIdx; idx < numCandidates; idx++)
	{
		/* Old Version (AB, AC, therefore ABC and ACB)
		if ( strncmp( (char*)&h_episodeCandidates[currentIdx*(prefixLength+1)],
					  (char*)&h_episodeCandidates[idx*(prefixLength+1)],
					  prefixLength ) == 0 && h_episodeSupport[idx] >= support )
		{
			return idx;
		}
		*/

		// New Version (AB, BC, therefore ABC)
		if ( strncmp( (char*)&h_episodeCandidates[currentIdx*(prefixLength+1)+1],
					  (char*)&h_episodeCandidates[idx*(prefixLength+1)],
					  prefixLength ) == 0 && h_episodeSupport[idx] >= support
					  && h_episodeCandidates[currentIdx*(prefixLength+1)+0] != h_episodeCandidates[idx*(prefixLength+1)+prefixLength])
		{
			bool intervalMatch = true;
			for( unsigned int intIdx = 0; intIdx < prefixLength-1; intIdx++ )
			{
				if ( h_episodeIntervals[currentIdx*prefixLength*2+(intIdx+1)*2+0] != h_episodeIntervals[idx*prefixLength*2+intIdx*2+0] ||
						 h_episodeIntervals[currentIdx*prefixLength*2+(intIdx+1)*2+1] != h_episodeIntervals[idx*prefixLength*2+intIdx*2+1] )
				{
					intervalMatch = false;
				}
			}

			if ( intervalMatch )
				return idx;
		}

	}

	return -1;
}

void triangleToArrayCPU( int triangle, int* base, int* compare )
{
	*base = 0;
	int Temp = triangle;
	while (Temp > 0)
	{
   		Temp = Temp - (numCandidates - *base);
   		(*base)++;
	}
	(*base)--;

	*compare = triangle + (*base-1) * (*base)/2 - (*base-1)*numCandidates;
}

void generateEpisodeCandidatesCPU( int level )
{
	int numCandidatesBuffer = 0;

	if ( level == 1 )
	{
		initEpisodeCandidates();
	}
	else if ( level == 2 )
	{
		for ( int candidateIdx1 = 0; candidateIdx1 < numCandidates; candidateIdx1++ )
		{
			if ( h_episodeSupport[candidateIdx1] < support )
					continue;

			for ( int candidateIdx2 = 0; candidateIdx2 < numCandidates; candidateIdx2++ )
			{
				if ( h_episodeSupport[candidateIdx2] < support || h_episodeCandidates[candidateIdx1] == h_episodeCandidates[candidateIdx2] )
					continue;

				for ( unsigned int idx = 0; idx < temporalConstraintSize; idx++ )
				{
					h_episodeCandidatesBuffer[numCandidatesBuffer*level+0] = h_episodeCandidates[candidateIdx1];
					h_episodeCandidatesBuffer[numCandidatesBuffer*level+1] = h_episodeCandidates[candidateIdx2];
					h_episodeIntervalsBuffer[numCandidatesBuffer*(level-1)*2+0] = temporalConstraint[idx*2+0];
					h_episodeIntervalsBuffer[numCandidatesBuffer*(level-1)*2+1] = temporalConstraint[idx*2+1];
					numCandidatesBuffer++;
				}
			}
		}


		UBYTE* tempCandidates = h_episodeCandidates;
		h_episodeCandidates = h_episodeCandidatesBuffer;
		h_episodeCandidatesBuffer = tempCandidates;
		numCandidates = numCandidatesBuffer;

		float* tempIntervals = h_episodeIntervals;
		h_episodeIntervals = h_episodeIntervalsBuffer;
		h_episodeIntervalsBuffer = tempIntervals;
	}
	else
	{
		for ( int candidateIdx = 0; candidateIdx < numCandidates-1; candidateIdx++ )
		{
			if ( numCandidates == 0 )
				break;

			if ( h_episodeSupport[candidateIdx] < support )
				continue;


			int nextCandidateIdx = 0;

			while ( (nextCandidateIdx = getNextValidCandidate(level-2, candidateIdx, nextCandidateIdx)) != -1)
			{
				strncpy( (char*)&h_episodeCandidatesBuffer[numCandidatesBuffer*level],
						 (char*)&h_episodeCandidates[candidateIdx*(level-1)],
						 level-1 );
				//h_episodeCandidatesBuffer[numCandidatesBuffer*level+level-2] =
				//	h_episodeCandidates[candidateIdx*(level-1)+level-2];
				h_episodeCandidatesBuffer[numCandidatesBuffer*level+level-1] =
					h_episodeCandidates[nextCandidateIdx*(level-1)+level-2];

				memcpy( &h_episodeIntervalsBuffer[numCandidatesBuffer*(level-1)*2],
						&h_episodeIntervals[candidateIdx*(level-2)*2],
						(level-1)*2*sizeof(float) );

				memcpy( &h_episodeIntervalsBuffer[numCandidatesBuffer*(level-1)*2+(level-2)*2],
								&h_episodeIntervals[nextCandidateIdx*(level-2)*2+(level-3)*2],
								2*sizeof(float) );
				//h_episodeIntervalsBuffer[numCandidatesBuffer*(level-1)*2+(level-2)*2+0] =
				//	h_episodeIntervals[candidateIdx*(level-2)*2+(level-3)*2+0];
				//h_episodeIntervalsBuffer[numCandidatesBuffer*(level-1)*2+(level-2)*2+1] =
				//	h_episodeIntervals[candidateIdx*(level-2)*2+(level-3)*2+1];

				numCandidatesBuffer++;
				nextCandidateIdx++;
			}
		}

		UBYTE* tempCandidates = h_episodeCandidates;
		h_episodeCandidates = h_episodeCandidatesBuffer;
		h_episodeCandidatesBuffer = tempCandidates;
		numCandidates = numCandidatesBuffer;

		float* tempIntervals = h_episodeIntervals;
		h_episodeIntervals = h_episodeIntervalsBuffer;
		h_episodeIntervalsBuffer = tempIntervals;
	}
}

void cullCandidates( int level )
{
	int numCandidatesBuffer = 0;
    printf("Culling Candidates\n");
	for ( int candidateIdx = 0; candidateIdx < numCandidates; candidateIdx++ )
	{
		if ( h_episodeSupport[candidateIdx] < support )
			continue;

		strncpy( (char*)&h_episodeCandidatesBuffer[numCandidatesBuffer*level],
				 (char*)&h_episodeCandidates[candidateIdx*level],
				 level );

		memcpy( &h_episodeIntervalsBuffer[numCandidatesBuffer*(level-1)*2],
				&h_episodeIntervals[candidateIdx*(level-1)*2],
				(level-1)*2*sizeof(float) );

		numCandidatesBuffer++;
	}

	episodesCulled = numCandidates - numCandidatesBuffer;

	UBYTE* tempCandidates = h_episodeCandidates;
	h_episodeCandidates = h_episodeCandidatesBuffer;
	h_episodeCandidatesBuffer = tempCandidates;
	numCandidates = numCandidatesBuffer;

	float* tempIntervals = h_episodeIntervals;
	h_episodeIntervals = h_episodeIntervalsBuffer;
	h_episodeIntervalsBuffer = tempIntervals;
}

void saveResult(int level)
{
	fprintf( dumpFile, "-----------------------\nEpisodes of size = %d\n-----------------------\n", level );

	unsigned int newcount = 0;
	char c; int v;
	for ( int idx = 0; idx < numCandidates; idx++ )
	{
		if ( h_episodeSupport[idx] >= support )
		{
			newcount++;
			for ( int levelIdx = 0; levelIdx < level; levelIdx++ )
			{
				if ( levelIdx > 0 )
				{
					fprintf(dumpFile, "-[%f,%f]-", h_episodeIntervals[idx*(level-1)*2+levelIdx*2+0], h_episodeIntervals[idx*(level-1)*2+levelIdx*2+1]);
				}
				if ( symbolSize == ONE_CHAR )
					fprintf( dumpFile, "%c", h_episodeCandidates[idx*level+levelIdx] );
				else
				{
					charToSymbol( h_episodeCandidates[idx*level+levelIdx], c, v );
					fprintf( dumpFile, "%c%d", c, v );
				}
			}
			fprintf( dumpFile, ": %f\n", h_episodeSupport[idx]);
		}
	}
	fprintf( dumpFile, "No. of %d node frequent episodes = %d\n\n", level, newcount);
}

unsigned int countLinesInFile( char* filename )
{
	int ch, prev = '\n' /* so empty files have no lines */, lines = 0;

	FILE* file = fopen( filename, "r" );

	// Count lines, each line is an symbol/timestamp pair
	// Retrieved from http://www.daniweb.com/code/snippet325.html
	// Author: Dave Sinkula
	if ( file )
	{
		while ( (ch = fgetc(file)) != EOF ) /* Read all chars in the file. */
		{
			if ( ch == '\n' )
			{
				++lines; /* Bump the counter for every newline. */
			}
			prev = ch; /* Keep a copy to later test whether... */
		}
		fclose(file);

		if ( prev != '\n' ) /* ...the last line did not end in a newline. */
		{
			++lines; /* If so, add one more to the total. */
		}
	}
	return lines;
}

int chooseAlgorithmType( int lev, long num, int threadsPerBlock )
{
	if ( lev > 10 )
		lev = 10;
	if ( lev == 1 )
		return MAP_AND_MERGE;
	else
		return NAIVE;

	if ( lev == 1)
		return MAP_AND_MERGE;
	else
		return NAIVE;
	int bpmp;
	const int maxBlocksPerMultiprocessor = 8;

	// Implementation as hypothesized by Dr. Cao
	const int mp = 30;
	if ( lev == 1 )
		bpmp = maxBlocksPerMultiprocessor;
	else
		bpmp = min(16384/((lev-1)*11*4),maxBlocksPerMultiprocessor);  // Blocks per multiprocessor is limited by
	                                               // shared mem per blockin this program, up to 8
	int tpb = threadsPerBlock;
  return num < mp*bpmp*tpb ? MAP_AND_MERGE : NAIVE;

	// Original Implementation
	//return num < optimalThreshold[lev-1] ? MAP_AND_MERGE : NAIVE; // Lookup table, from OptimalThreshold.h
}

int loadData(char* filename)
{
	FILE* eventFile;

	eventSize = countLinesInFile(filename);
	if ((eventFile =  fopen( filename, "r" )) == NULL )
	{
		printf("Can not open data file: %s\n", filename);
		return -1;
	}

	int remainder = eventSize % MaxSections;
	int toPad = MaxSections - remainder;

	padEventSize = eventSize + toPad;

	h_events = (UBYTE*)malloc(padEventSize * sizeof(UBYTE));
	h_times = (float*)malloc(padEventSize * sizeof(float));

	// test file for one or two-char inputs
	char c1, c2;
	fscanf( eventFile, "%c%c", &c1, &c2 );
	if ( c2 == ',' )
		symbolSize = ONE_CHAR;
	else
		symbolSize = TWO_CHAR;
	rewind( eventFile );

	char symbol = 0;
	char c; int v;
	float time = 0;
	for ( unsigned int idx = 0; idx < eventSize; idx++ )
	{
		if ( symbolSize == ONE_CHAR )
			fscanf( eventFile, "%c,%f\n", &symbol, &time );
		else
		{
			fscanf( eventFile, "%c%d,%f\n", &c, &v, &time );
			symbol = symbolToChar( c, v );
		}
		h_events[idx] = symbol;
		h_times[idx] = time;
	}

	memset(&h_events[eventSize], 0, toPad*sizeof(UBYTE));

	for ( unsigned int idx = eventSize; idx < padEventSize; idx++ )
	{
		h_times[idx] = 1000000000.0f;
	}
	fclose( eventFile );
	return 0;
}

int loadTemporalConstraints(char* filename)
{
	FILE* temporalFile;

	temporalConstraintSize = countLinesInFile(filename);
	if ( (temporalFile = fopen( filename, "r" )) == NULL)
	{
		printf("Can not open temporal constraints file: %s\n", filename);
		return -1;
	}

	temporalConstraint = (float*)malloc(temporalConstraintSize*2);

	for ( unsigned int idx = 0; idx < temporalConstraintSize; idx++ )
	{
		fscanf( temporalFile, "%f %f\n", &temporalConstraint[2*idx+0], &temporalConstraint[2*idx+1] );
	}
	return 0;
}

void initGpu()
{
    int err;
    /////////////////////////////////////////////////////////////
    // Basic OpenCL Setup

	device_id = GetDevice(platform_id, n_device);

    // Create a compute context
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");

    // Create a command queue
    commands = clCreateCommandQueue(context, device_id, TIMER_ENABLE, &err);
    CHKERR(err, "Failed to create a command queue!");
    /////////////////////////////////////////////////////////////
}

void setupGpu()
{
    int err;


    /////////////////////////////////////////////////////////////
    //Compile Source

    //Read file
    FILE* kernelFile = fopen("GpuTemporalDataMining.cl", "r");
    fseek(kernelFile, 0, SEEK_END);
    size_t kernelLength = (size_t) ftell(kernelFile);
    char* kernelSource = (char *) malloc(sizeof(char)*(kernelLength+1));
    rewind(kernelFile);
    int read = fread((void *) kernelSource, kernelLength, 1, kernelFile);
    fclose(kernelFile);
    kernelSource[kernelLength] = 0;
     // Create the compute program from the source buffer
    program = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, NULL, &err);
    CHKERR(err, "Failed to create a compute program!");
    
	printf("MaxImageWidth: %d, MaxImageHeight: %d\n",(int)MaxImageWidth, (int)MaxImageHeight);
    // Build the program executable
    char options[200];
    snprintf(options, 200, "-I . -D IMAGE_MAX_WIDTH=%lu -D IMAGE_MAX_HEIGHT=%lu", (unsigned long)MaxImageWidth, (unsigned long)MaxImageHeight);
    err = clBuildProgram(program, 0, NULL, options, NULL, NULL);
    if (err == CL_BUILD_PROGRAM_FAILURE)
    {
        char *log;
        size_t logLen;
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLen);
        log = (char *) malloc(sizeof(char)*logLen);
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logLen, (void *) log, NULL);
        fprintf(stderr, "CL Error %d: Failed to build program! Log:\n%s", err, log);
        free(log);
        exit(1);
    }
    CHKERR(err, "Failed to build program!");

    free(kernelSource);

    kernel_countCandidates = clCreateKernel(program, "countCandidates", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    kernel_countCandidatesStatic = clCreateKernel(program, "countCandidatesStatic", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    kernel_countCandidatesMapMerge = clCreateKernel(program, "countCandidatesMapMerge", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    kernel_countCandidatesMapMergeStatic = clCreateKernel(program, "countCandidatesMapMergeStatic", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    kernel_analyzeSupport = clCreateKernel(program, "analyzeSupport", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    kernel_generateEpisodeCandidatesKernel = clCreateKernel(program, "generateEpisodeCandidatesKernel", &err);
    CHKERR(err, "Failed to create a compute kernel!");
    /////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////
    //Setup memory


    d_events = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, roundup(padEventSize) * sizeof(UBYTE), (void*)h_events, &err);
    CHKERR(err, "Failed to allocate device memory!");
    bindTexture(0, &eventTex, d_events, roundup(padEventSize) * sizeof(UBYTE), CL_UNSIGNED_INT8);

	free( h_events );

	d_times = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, roundup(padEventSize) * sizeof(float), (void*)h_times, &err);
    CHKERR(err, "Failed to allocate device memory!");
    bindTexture(0, &timeTex, d_times, roundup(padEventSize), CL_FLOAT);
	free( h_times );


    d_episodeCandidates = clCreateBuffer(context, CL_MEM_READ_WRITE, maxCandidates * sizeof(UBYTE), NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    d_episodeCandidatesBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, maxCandidates * sizeof(UBYTE), NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");

    h_episodeCandidates = (UBYTE*)malloc( maxCandidates*sizeof(UBYTE) );
	h_episodeCandidatesBuffer = (UBYTE*)malloc( maxCandidates*sizeof(UBYTE) );


	d_episodeIntervals = clCreateBuffer(context, CL_MEM_READ_WRITE, maxIntervals * sizeof(float), NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
	d_episodeIntervalsBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, maxIntervals * sizeof(float), NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");

	h_episodeIntervals = (float*)malloc( maxIntervals*sizeof(float) );
	h_episodeIntervalsBuffer = (float*)malloc( maxIntervals*sizeof(float) );

    // Results
	h_episodeSupport = (float*)malloc( maxCandidates*sizeof(float) );
	d_episodeSupport = clCreateBuffer(context, CL_MEM_READ_WRITE, maxCandidates * sizeof(float), NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");

	//h_mapRecords = (float*)malloc( 3 * numSections * maxLevel * maxCandidates * sizeof(float) );
	//CUDA_SAFE_CALL( cudaMalloc( (void**)&d_mapRecords, 3 * numSections * maxLevel * maxCandidates * sizeof(float)) );

    /////////////////////////////////////////////////////////////
}

void countCandidates(size_t* globalWork, size_t* localWork, cl_mem episodeSupport, long eventSize, int level, int sType, int numCandidates,
                    cl_mem candidateTex, cl_mem intervalTex, cl_mem eventTex, cl_mem timeTex, size_t sharedMemNeeded)
{
    int errcode;
    errcode  = clSetKernelArg(kernel_countCandidates, 0, sizeof(cl_mem), (void *) &episodeSupport);
    errcode |= clSetKernelArg(kernel_countCandidates, 1, sizeof(long), (void *) &eventSize);
    errcode |= clSetKernelArg(kernel_countCandidates, 2, sizeof(int), (void *) &level);
    errcode |= clSetKernelArg(kernel_countCandidates, 3, sizeof(int), (void *) &sType);
    errcode |= clSetKernelArg(kernel_countCandidates, 4, sizeof(int), (void *) &numCandidates);
    errcode |= clSetKernelArg(kernel_countCandidates, 5, sizeof(cl_mem), (void *) &candidateTex);
    errcode |= clSetKernelArg(kernel_countCandidates, 6, sizeof(cl_mem), (void *) &intervalTex);
    errcode |= clSetKernelArg(kernel_countCandidates, 7, sizeof(cl_mem), (void *) &eventTex);
    errcode |= clSetKernelArg(kernel_countCandidates, 8, sizeof(cl_mem), (void *) &timeTex);
    errcode |= clSetKernelArg(kernel_countCandidates, 9, sharedMemNeeded, NULL);
    CHKERR(errcode, "Unable to set arguments for countCandidates");
    START_TIMER
	errcode = clEnqueueNDRangeKernel(commands, kernel_countCandidates, 3, NULL, globalWork, localWork, 0, NULL, &myEvent);
    	clFinish(commands);
	END_TIMER
	COUNT_K
    CHKERR(errcode, "Error running countCandidates");
}

void countCandidatesMapMerge(size_t* globalWork, size_t* localWork, cl_mem episodeSupport, long padEventSize, int level, int sType, int sections, int x, int numCandidates,
                    cl_mem candidateTex, cl_mem intervalTex, cl_mem eventTex, cl_mem timeTex, size_t sharedMemNeeded)
{
    int errcode;
    errcode  = clSetKernelArg(kernel_countCandidatesMapMerge, 0, sizeof(cl_mem), (void *) &episodeSupport);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 1, sizeof(long), (void *) &padEventSize);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 2, sizeof(int), (void *) &level);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 3, sizeof(int), (void *) &sType);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 4, sizeof(int), (void *) &sections);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 5, sizeof(int), (void *) &x);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 6, sizeof(int), (void *) &numCandidates);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 7, sizeof(cl_mem), (void *) &candidateTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 8, sizeof(cl_mem), (void *) &intervalTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 9, sizeof(cl_mem), (void *) &eventTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 10,sizeof(cl_mem), (void *) &timeTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMerge, 11, sharedMemNeeded, NULL);
    CHKERR(errcode, "Unable to set arguments for countCandidatesMapMerge");
    START_TIMER
	errcode = clEnqueueNDRangeKernel(commands, kernel_countCandidatesMapMerge, 3, NULL, globalWork, localWork, 0, NULL, &myEvent);
    	clFinish(commands);
	END_TIMER
	COUNT_K
	CHKERR(errcode, "Error running countCandidatesMapMerge");
}
void countCandidatesStatic(size_t* globalWork, size_t* localWork, cl_mem episodeSupport, long eventSize, int level, int sType, int numCandidates,
                    cl_mem candidateTex, cl_mem intervalTex, cl_mem eventTex, cl_mem timeTex, size_t sharedMemNeeded)
{
    int errcode;
    errcode  = clSetKernelArg(kernel_countCandidatesStatic, 0, sizeof(cl_mem), (void *) &episodeSupport);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 1, sizeof(long), (void *) &eventSize);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 2, sizeof(int), (void *) &level);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 3, sizeof(int), (void *) &sType);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 4, sizeof(int), (void *) &numCandidates);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 5, sizeof(cl_mem), (void *) &candidateTex);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 6, sizeof(cl_mem), (void *) &intervalTex);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 7, sizeof(cl_mem), (void *) &eventTex);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 8, sizeof(cl_mem), (void *) &timeTex);
    errcode |= clSetKernelArg(kernel_countCandidatesStatic, 9, sharedMemNeeded, NULL);
    CHKERR(errcode, "Unable to set arguments for countCandidates");
    START_TIMER
	errcode = clEnqueueNDRangeKernel(commands, kernel_countCandidatesStatic, 3, NULL, globalWork, localWork, 0, NULL, &myEvent);
    	clFinish(commands);
	END_TIMER
	COUNT_K
    CHKERR(errcode, "Error running countCandidatesStatic");
}

void countCandidatesMapMergeStatic(size_t* globalWork, size_t* localWork, cl_mem episodeSupport, long padEventSize, int level, int sType, int sections, int x, int numCandidates,
                    cl_mem candidateTex, cl_mem intervalTex, cl_mem eventTex, cl_mem timeTex, size_t sharedMemNeeded)
{
    int errcode;
    errcode  = clSetKernelArg(kernel_countCandidatesMapMergeStatic, 0, sizeof(cl_mem), (void *) &episodeSupport);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 1, sizeof(long), (void *) &padEventSize);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 2, sizeof(int), (void *) &level);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 3, sizeof(int), (void *) &sType);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 4, sizeof(int), (void *) &sections);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 5, sizeof(int), (void *) &x);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 6, sizeof(int), (void *) &numCandidates);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 7, sizeof(cl_mem), (void *) &candidateTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 8, sizeof(cl_mem), (void *) &intervalTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 9, sizeof(cl_mem), (void *) &eventTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 10,sizeof(cl_mem), (void *) &timeTex);
    errcode |= clSetKernelArg(kernel_countCandidatesMapMergeStatic, 11, sharedMemNeeded, NULL);
    CHKERR(errcode, "Unable to set arguments for countCandidatesMapMerge");
    START_TIMER
	errcode = clEnqueueNDRangeKernel(commands, kernel_countCandidatesMapMergeStatic, 3, NULL, globalWork, localWork, 0, NULL, &myEvent);
    	clFinish(commands);
	END_TIMER
	COUNT_K
    CHKERR(errcode, "Error running countCandidatesMapMergeStatic");
}

void cleanup()
{
	printf("Cleaning up memory...\n");
    clReleaseMemObject(eventTex);
    clReleaseMemObject(timeTex);

	clReleaseMemObject(d_events);
	clReleaseMemObject(d_episodeSupport);
	clReleaseMemObject(d_episodeCandidates);
	clReleaseMemObject(d_episodeCandidatesBuffer);

	free( h_episodeSupport );
	free( h_episodeCandidates );
	free( h_episodeCandidatesBuffer );

	fclose(dumpFile);
}

void calculateGrid(size_t grid[3], int threads, int candidates )
{
	//dim3 grid;
	////if ( numCandidates / threads <= 16 )
	////{
	////	int temp = numCandidates / 16 + 1;
	////	grid = dim3(numCandidates / temp + 1, 1, 1);
	////}
	////else
	////{
	////	grid = dim3( numCandidates / threads + 1, 1, 1);
	////}
	//grid = dim3( numCandidates / threads + (numCandidates % threads == 0 ? 0 : 1), 1, 1);
	//return grid;
	grid[0] = threads * (numCandidates / threads + (numCandidates % threads == 0 ? 0 : 1));
    //printf("t = %d, numCan = %d\n", threads, numCandidates);
    //printf("g0 = %d\n", grid[0]);
	grid[1] = 1;
	grid[2] = 1;
}

void calculateBlock(size_t block[3], int threads, int candidates )
{
	//dim3 block;
	////if ( numCandidates / threads <= 16 )
	////{
	////	int temp = numCandidates / 16 + 1;
	////	block = dim3( temp, 1, 1);
	////}
	////else
	////{
	////	block = dim3( threads, 1, 1);
	////}

	//block = dim3( threads, 1, 1);
	//return block;
	block[0] = threads;
	block[1] = 1;
	block[2] = 1;
}

void getDeviceVariables(cl_device_id device)
{
	MinThreads = 32;
	int err;
    err = clGetDeviceInfo(device_id,CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&MaxThreads, NULL);
    CHKERR(err, "Error checking for work group size\n");
	err = clGetDeviceInfo(device_id,CL_DEVICE_LOCAL_MEM_SIZE,sizeof(size_t),&MaxSharedMemory, NULL);
    CHKERR(err, "Error checking for local memory size\n");

	size_t maxThreads[3];
	err = clGetDeviceInfo(device_id,CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t)*3,&maxThreads, NULL);
    CHKERR(err, "Error checking for work item sizes\n");
	MaxSections = maxThreads[0];

    cl_bool images;
	err = clGetDeviceInfo(device_id,CL_DEVICE_IMAGE_SUPPORT,sizeof(cl_bool),&images, NULL);
    CHKERR(err, "Error checking for image capability\n");
    if(images == CL_FALSE)
    {
        printf("This device does not support images\n");
        exit(1);
    }
	err = clGetDeviceInfo(device_id,CL_DEVICE_IMAGE2D_MAX_WIDTH,sizeof(size_t),&MaxImageWidth, NULL);
    CHKERR(err, "Unable to get image width\n");

	err = clGetDeviceInfo(device_id,CL_DEVICE_IMAGE2D_MAX_HEIGHT,sizeof(size_t),&MaxImageHeight, NULL);
    CHKERR(err, "Unable to get image width\n");

    //printf("MaxImageWidth: %d\n", MaxImageWidth);
    //printf("MaxImageHeight: %d\n", MaxImageHeight);

    //printf( "MinThreads: %d\n", MinThreads );
	//printf( "MaxThreads: %d\n", MaxThreads );
	//printf( "MaxSharedMemory: %d\n", MaxSharedMemory );
	//printf( "MaxSections: %d\n", MaxSections );
}

void calculateLevelParameters(int level, size_t* block, size_t* grid, int& sections)
{
	const int maxDim = 65536;
	MaxSections = 32;
	sections = (int)floor( (float)(MaxSharedMemory-384) / (float)(4*(level-1)*level*MaxListSize+12*level) );
	sections = sections > MaxSections ? MaxSections : sections;

	sections = (int)pow(2,floor(log((double)sections)/log(2.0)));

	int y = (numCandidates / maxDim);
	int lol = numCandidates % maxDim == 0 ? 0 : 1;
	grid[0] = sections * (numCandidates > maxDim ? maxDim : numCandidates);
	grid[1] = level * (y + lol);
	grid[2] = 1;
	//grid = dim3( numCandidates > maxDim ? maxDim : numCandidates,
	//			y + lol,
	//			1 );
	block[0] = sections;
	block[1] = level;
	block[2] = 1;
	//block = dim3( sections, level, 1 );
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv)
{
    INI_TIMER
	runTest( argc, argv);
	PRINT_COUNT
}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void
runTest( int argc, char** argv)
{
	if ( argc != 8 && argc != 10)
	{
		printf("Usage: GpuTemporalDataMining <file path> <temporal constraint path> <threads> <support> <(a)bsolute or (r)atio> <(s)tatic or (d)ynamic> <(m)ap and merge or (n)aive or (o)hybrid> [platform & device]\n");
		return;
	}

    if(argc == 10)
	{
		platform_id = atoi(argv[8]);
		n_device = atoi(argv[9]);
	}
//    CUT_DEVICE_INIT();
    initGpu();

    getDeviceVariables(device_id);


	printf("Dataset, Support Threshold, PTPE or MapMerge, A1 or A1+A2, Level, Episodes (N), Episodes Culled (X), A1 Counting Time, A2 Counting Time, Generation Time, Total Counting Time\n");

    //CUT_SAFE_CALL( cutCreateTimer( &timer));
    //CUT_SAFE_CALL( cutCreateTimer( &generating_timer));
    //CUT_SAFE_CALL( cutCreateTimer( &a1_counting_timer));
    //CUT_SAFE_CALL( cutCreateTimer( &a2_counting_timer));
    //CUT_SAFE_CALL( cutCreateTimer( &total_timer));

    //CUT_SAFE_CALL( cutStartTimer( total_timer));
    //CUT_SAFE_CALL( cutStartTimer( timer));
    //CUT_SAFE_CALL( cutStartTimer( generating_timer));
    //CUT_SAFE_CALL( cutStartTimer( a1_counting_timer));
    //CUT_SAFE_CALL( cutStartTimer( a2_counting_timer));

    unsigned int num_threads = atoi(argv[3]);
    // allocate host memory
    //initEpisodeCandidates();
	if ( loadData( argv[1] ) != 0 )
		return;
	if ( loadTemporalConstraints(argv[2]) != 0 )
		return;

	// Check whether value supplied is absolute or ratio support
	supportType = *(argv[5]) == 'a' ? ABSOLUTE : RATIO;
	memoryModel = *(argv[6]) == 's' ? STATIC : DYNAMIC;

	switch (*(argv[7]))
	{
	case 'm':
		algorithmType = MAP_AND_MERGE;
		break;
	case 'n':
		algorithmType = NAIVE;
		break;
	case 'o':
		algorithmType = OPTIMAL;
		break;
	}

	support = atof(argv[4]);

	dumpFile = fopen( "episode.txt", "w" );

	//printf("Initializing GPU Data...\n");

	setupGpu();

	// setup execution parameters
    size_t grid[3];
    size_t threads[3];

	//printf("Event stream size: %i\n", eventSize);

	// BEGIN LOOP
	for ( int level = 1; level <= eventSize; level++ )
	{
		printf("Generating episode candidates for level %i...\n", level);
//		CUT_SAFE_CALL( cutResetTimer( total_timer));
//		CUT_SAFE_CALL( cutStartTimer( total_timer));

		//CUDA_SAFE_CALL( cudaUnbindTexture( candidateTex ) );
		if(level != 1){
			unbindTexture(&candidateTex, d_episodeCandidates, numCandidates * (level-1) * sizeof(UBYTE) );
		//CUDA_SAFE_CALL( cudaUnbindTexture( intervalTex ) );
		unbindTexture(&intervalTex, d_episodeIntervals, numCandidates * (level-2) * 2 * sizeof(float));
        }

//		CUT_SAFE_CALL( cutResetTimer( generating_timer));
//		CUT_SAFE_CALL( cutStartTimer( generating_timer));

//		int test1, test = numCandidates;
//		generateEpisodeCandidatesCPU( level );
//		test1 = numCandidates;
//		numCandidates = test;
        printf("Generating Episodes\n");
#ifdef CPU_EPISODE_GENERATION
		generateEpisodeCandidatesCPU( level );
#else
		generateEpisodeCandidatesGPU( level, num_threads );
#endif

//		CUT_SAFE_CALL( cutStopTimer( generating_timer));
		//printf( "\tGenerating time: %f (ms)\n", cutGetTimerValue( generating_timer));


		if ( numCandidates == 0 )
			break;
        printf("Writing to buffer\n");
		// Copy candidates to GPU
#ifdef CPU_EPISODE_GENERATION
		START_TIMER
		clEnqueueWriteBuffer(commands, d_episodeCandidates, CL_TRUE, 0, numCandidates * level * sizeof(UBYTE), h_episodeCandidates, 0, NULL, &myEvent);
		END_TIMER
		COUNT_H2D
		clEnqueueWriteBuffer(commands, d_episodeIntervals, CL_TRUE, 0, numCandidates * (level-1) * 2 * sizeof(float), h_episodeIntervals, 0, NULL, &myEvent);
	CL_FINISH(commands)
		END_TIMER
		COUNT_H2D
#endif

        bindTexture( 0, &candidateTex, d_episodeCandidates, numCandidates * level * sizeof(UBYTE), CL_UNSIGNED_INT8);
        
		bindTexture( 0, &intervalTex, d_episodeIntervals, numCandidates * (level-1) * 2 * sizeof(float), CL_FLOAT );

		//printf("Executing kernel on %i candidates...\n", numCandidates, level);

		// execute the kernel
		calculateGrid(grid, num_threads, numCandidates);
		calculateBlock(threads, num_threads, numCandidates);

		int sections;
		unsigned int shared_mem_needed;
		//CUT_SAFE_CALL( cutStartTimer( counting_timer));

		int aType = algorithmType;
		if ( algorithmType == OPTIMAL )
			aType = chooseAlgorithmType( level, numCandidates, num_threads );

		if ( memoryModel == DYNAMIC )
		{
			if ( aType == NAIVE )
			{
				shared_mem_needed = MaxListSize*level*threads[0]*sizeof(float);
                printf("Shared memory needed %d\n", shared_mem_needed);
				//CUT_SAFE_CALL( cutResetTimer( a1_counting_timer));
				//CUT_SAFE_CALL( cutStartTimer( a1_counting_timer));
				countCandidates(grid, threads, d_episodeSupport, eventSize, level, supportType, numCandidates, candidateTex, intervalTex, eventTex, timeTex, shared_mem_needed );

			}
			else
			{
                printf("DYNAMIC MAP MERGE\n");
				calculateLevelParameters(level, threads, grid, sections);
				shared_mem_needed = 16000;
                printf("numCandidates=%d\n", numCandidates);
				//CUT_SAFE_CALL( cutResetTimer( a1_counting_timer));
				//CUT_SAFE_CALL( cutStartTimer( a1_counting_timer));
				countCandidatesMapMerge(grid, threads, d_episodeSupport, padEventSize, level, supportType, sections, padEventSize / sections, numCandidates,
                    candidateTex, intervalTex, eventTex, timeTex, shared_mem_needed );
				//countCandidatesMapMergeStatic<<< grid, threads, shared_mem_needed >>>( d_episodeSupport, padEventSize, level, supportType, sections, padEventSize / sections, numCandidates );
			}
		}
		else
		{
			if ( aType == NAIVE )
			{
				shared_mem_needed = level*threads[0]*sizeof(float);
			}
			else
			{
				calculateLevelParameters(level, threads, grid, sections);
				shared_mem_needed = 16000;
			}

				//CUT_SAFE_CALL( cutResetTimer( a2_counting_timer));
				//CUT_SAFE_CALL( cutStartTimer( a2_counting_timer));
			if ( aType == NAIVE )
                countCandidatesStatic(grid, threads, d_episodeSupport, eventSize, level, supportType, numCandidates, candidateTex, intervalTex, eventTex, timeTex, shared_mem_needed  );
			else
                countCandidatesMapMergeStatic(grid, threads, d_episodeSupport, padEventSize, level, supportType, sections, padEventSize / sections, numCandidates, 
                    candidateTex, intervalTex, eventTex, timeTex, shared_mem_needed );
			clFinish(commands);

			//CUT_SAFE_CALL( cutStopTimer( a2_counting_timer));

            int err;
		START_TIMER
            err = clEnqueueReadBuffer(commands,d_episodeSupport, CL_TRUE, 0, numCandidates * sizeof(float), h_episodeSupport, 0, NULL, &myEvent);
	CL_FINISH(commands)
            	END_TIMER
		COUNT_D2H
		CHKERR(err, "Unable to read buffer from device.");
		
			unbindTexture(&candidateTex, d_episodeCandidates, numCandidates * level * sizeof(UBYTE) );
			unbindTexture(&intervalTex, d_episodeIntervals, numCandidates * (level-1) * 2 * sizeof(float));

			// Remove undersupported episodes
			cullCandidates( level );

			if ( numCandidates == 0 )
				break;
			unsigned int mmthreads = num_threads;
			if ( MaxListSize*level*num_threads*sizeof(float) > 16384 )
			{
				if ( MaxListSize*level*96*sizeof(float) < 16384 )
					mmthreads = 96;
				else if ( MaxListSize*level*64*sizeof(float) < 16384)
					mmthreads = 64;
				else if ( MaxListSize*level*32*sizeof(float) < 16384)
					mmthreads = 32;
				printf("More shared memory needed for %d threads. Changed to %d threads.\n", num_threads, mmthreads );
			}

#ifdef CPU_EPISODE_GENERATION
		START_TIMER
            err = clEnqueueWriteBuffer(commands, d_episodeCandidates, CL_TRUE, 0, numCandidates * level * sizeof(UBYTE), h_episodeCandidates, 0, NULL, &myEvent);
		END_TIMER
		COUNT_H2D
            CHKERR(err, "Unable to write buffer 1.");
            if(numCandidates * (level - 1) * 2 * sizeof(float) != 0)
            err = clEnqueueWriteBuffer(commands, d_episodeIntervals, CL_TRUE, 0, numCandidates * (level-1) * 2 * sizeof(float), h_episodeIntervals, 0, NULL, &myEvent);
	CL_FINISH(commands)
            CHKERR(err, "Unable to write buffer 2.");
		END_TIMER
		COUNT_H2D
#endif
            bindTexture( 0, &candidateTex, d_episodeCandidates, numCandidates * level * sizeof(UBYTE), CL_UNSIGNED_INT8);
            bindTexture( 0, &intervalTex, d_episodeIntervals, numCandidates * (level-1) * 2 * sizeof(float), CL_FLOAT );

			if ( algorithmType == OPTIMAL )
				aType = chooseAlgorithmType( level, numCandidates, mmthreads );

			// Run (T1,T2] algorithm
			if ( aType == NAIVE )
			{
				shared_mem_needed = MaxListSize*level* mmthreads*sizeof(float);
				calculateGrid(grid, mmthreads, numCandidates );
				calculateBlock(threads, mmthreads, numCandidates );
			}
			else
			{
				calculateLevelParameters(level, threads, grid, sections);
				shared_mem_needed = 16000;
			}
			//CUT_SAFE_CALL( cutResetTimer( a1_counting_timer));
			//CUT_SAFE_CALL( cutStartTimer( a1_counting_timer));

			if ( aType == NAIVE )
                countCandidates(grid, threads, d_episodeSupport, eventSize, level, supportType, numCandidates,
                    candidateTex, intervalTex, eventTex, timeTex, shared_mem_needed );
			else
                countCandidatesMapMerge(grid, threads, d_episodeSupport, padEventSize, level, supportType, sections, padEventSize / sections, numCandidates,
                    candidateTex, intervalTex, eventTex, timeTex, shared_mem_needed );

		}
        printf("Finishing\n");
		clFinish(commands);
		//CUT_SAFE_CALL( cutStopTimer( a1_counting_timer));
		//printf( "\tCounting time: %f (ms)\n", cutGetTimerValue( counting_timer));

		// check if kernel execution generated an error
		//CUT_CHECK_ERROR("Kernel execution failed");

		//printf("Copying result back to host...\n\n");
		START_TIMER
        int err = clEnqueueReadBuffer(commands, d_episodeSupport, CL_TRUE, 0, numCandidates * sizeof(float), h_episodeSupport, 0, NULL, &myEvent);
	CL_FINISH(commands)
		END_TIMER
		COUNT_D2H
		CHKERR(err, "Unable to read memory 1.");
		err = clEnqueueReadBuffer(commands, d_episodeCandidates, CL_TRUE, 0, numCandidates * level * sizeof(UBYTE), h_episodeCandidates, 0, NULL, &myEvent);
	CL_FINISH(commands)
		CHKERR(err, "Unable to read memory 2.");
		//CUDA_SAFE_CALL( cudaMemcpy( h_mapRecords, d_mapRecords, 3 * numSections * maxLevel * maxCandidates * sizeof(float), cudaMemcpyDeviceToHost ));
		END_TIMER
		COUNT_D2H
		saveResult(level);
		fflush(dumpFile);
	// END LOOP

		//CUT_SAFE_CALL( cutStopTimer( total_timer));

		// Print Statistics for this run
		printf("%s, %f, %s, %s, %d, %d, %d\n",
			argv[1],											// Dataset
			support,											// Support Threshold
			algorithmType == NAIVE ? "PTPE" : algorithmType == MAP_AND_MERGE ? "MapMerge" : "Episode-Based",		// PTPE or MapMerge or Episode-Based
			memoryModel == STATIC ? "A1+A2" : "A1",				// A1 or A1+A2
			level,												// Level
			numCandidates+episodesCulled,						// Episodes counted
			episodesCulled  									// Episodes removed by A2
	//		cutGetTimerValue( a1_counting_timer),				// Time for A1
//			memoryModel == STATIC ? cutGetTimerValue( a2_counting_timer) : 0.0f,				// Time for A2
		//	cutGetTimerValue( generating_timer),				// Episode generation time
		//	cutGetTimerValue( total_timer) );					// Time for total loop
            );
	}
	printf("Done!\n");

    cleanup();

    //CUT_SAFE_CALL( cutStopTimer( timer));
    //printf( "Processing time: %f (ms)\n", cutGetTimerValue( timer));
    //CUT_SAFE_CALL( cutDeleteTimer( timer));
    //CUT_SAFE_CALL( cutDeleteTimer( generating_timer));
    //CUT_SAFE_CALL( cutDeleteTimer( a1_counting_timer));
    //CUT_SAFE_CALL( cutDeleteTimer( a2_counting_timer));
    //CUT_SAFE_CALL( cutDeleteTimer( total_timer));

}
