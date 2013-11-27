#ifndef GLOBAL_H_
#define GLOBAL_H_

const unsigned int MaxRecords = 5000000;
const unsigned int maxCandidates = 13000000;
//const unsigned int maxLevel = 20;
const unsigned int maxIntervals = (maxCandidates-1)*2;
const unsigned int MaxListSize = 10;

static float support;

enum
{
	EVENT_26,
	EVENT_64,
};

// Device Dependent Values
static int MaxSharedMemory, MinThreads, MaxThreads, MaxSections;

#endif /* GLOBAL_H_ */
