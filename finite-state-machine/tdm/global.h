#ifndef GLOBAL_H_
#define GLOBAL_H_

//MACRO OPTIONS

#define CPU_EPISODE_GENERATION

//const unsigned int maxWidth = 4096;
//const unsigned int maxHeight = 4096;
const unsigned int maxCandidates = 130000;
const unsigned int maxLevel = 20;
const unsigned int maxIntervals = (maxCandidates-1)*2;
const unsigned int MaxListSize = 11;
const unsigned int uniqueEvents = 'Z'-'A'+1;

//texture<UBYTE, 1, cudaReadModeElementType> candidateTex;
//texture<float, 1, cudaReadModeElementType> intervalTex;

cl_mem candidateTex;
cl_mem intervalTex;

//texture<UBYTE, 1, cudaReadModeElementType> eventTex;
//texture<float, 1, cudaReadModeElementType> timeTex;

cl_mem eventTex;
cl_mem timeTex;

static float support;

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

static int supportType, memoryModel, algorithmType;

// Device Dependent Values
static int MaxSharedMemory, MinThreads, MaxThreads, MaxSections;

#endif /* GLOBAL_H_ */
