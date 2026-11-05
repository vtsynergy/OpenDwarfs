#ifndef DATAIO_H_
#define DATAIO_H_

#include "types.h"
#include <stdio.h>

void loadData(char* filename, ubyte** events, float** times, uint* eventCount, uint* eventType, uint* uniqueEvents);
void loadTemporalConstraints(char* filename, float** constraints, uint* constraintCount);
void saveResult(FILE* dumpFile, int level, uint episodes, uint* support, ubyte* candidates, float* intervals, uint eventType);
char symbolToChar( char l, int n );
void charToSymbol( char c, char& l, int& n );

#endif /* DATAIO_H */
