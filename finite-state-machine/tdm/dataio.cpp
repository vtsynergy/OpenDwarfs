#include "dataio.h"
#include "global.h"

#include <stdio.h>
#include <stdlib.h>

char symbolToChar( char l, int n )
{
	return (l - 'a')*8 + (n-1) + '!';
}

void charToSymbol( char c, char& l, int& n )
{
	n = (c-'!')%8+1;
	l = (c-'!'-(n-1))/8 + 'a';	
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


void loadTemporalConstraints(char* filename, float** constraints, unsigned int* constraintCount)
{
	FILE* temporalFile;

	uint count = countLinesInFile(filename);
	temporalFile = fopen( filename, "r" );
	if ( !temporalFile ) {
	  printf("ERROR: %s does not exist.\n", filename);
	  return;
	}

	float* tc = (float*)malloc(count*2*sizeof(float));
	//printf("Count: %d\n", count);
	for ( uint idx = 0; idx < count; idx++ ) {
	  fscanf( temporalFile, "%f %f\n", &tc[2*idx+0], &tc[2*idx+1] );
	}
	fclose(temporalFile);

	*constraints = tc;
	*constraintCount = count;
}


void loadData(char* filename, ubyte** events, float** times, uint* eventCount, uint* eventType, uint* uniqueEvents)
{
	FILE* eventFile;

	uint eventSize;
	*eventCount = eventSize = countLinesInFile(filename);
	eventFile =  fopen( filename, "r" );
	
	*events = (ubyte*)malloc(eventSize * sizeof(ubyte));
	*times = (float*)malloc(eventSize * sizeof(float));

	// test file for one or two-char inputs
	char c1, c2;
	fscanf( eventFile, "%c%c", &c1, &c2 );
	if ( c2 == ',' ) {
		*eventType = EVENT_26;
		*uniqueEvents = 26;
	}
	else {
		*eventType = EVENT_64;
		*uniqueEvents = 64;
	}
	rewind( eventFile );
	
	char symbol = 0;
	char c; int v;
	float time = 0;
	for ( uint idx = 0; idx < eventSize; idx++ )
	{
		if ( *eventType == EVENT_26 )
			fscanf( eventFile, "%c,%f\n", &symbol, &time );
		else
		{
			fscanf( eventFile, "%c%d,%f\n", &c, &v, &time );
			symbol = symbolToChar( c, v );
		}
		(*events)[idx] = symbol;
		(*times)[idx] = time;
	}

	fclose( eventFile );
}


void saveResult(FILE* dumpFile, int level, uint episodes, uint* support, ubyte* candidates, float* intervals, uint eventType)
{
	fprintf( dumpFile, "-----------------------\nEpisodes of size = %d\n-----------------------\n", level );

	unsigned int newcount = 0;
	char c; int v;
	for ( int idx = 0; idx < episodes; idx++ )
	{
		newcount++;
		for ( int levelIdx = 0; levelIdx < level; levelIdx++ )
		{
			if ( levelIdx > 0 )
			{
				fprintf(dumpFile, "-[%f,%f]-", intervals[idx*(level-1)*2+levelIdx*2+0], intervals[idx*(level-1)*2+levelIdx*2+1]);
			}
			if ( eventType == EVENT_26 )
				fprintf( dumpFile, "%c", candidates[idx*level+levelIdx] );
			else
			{
				charToSymbol( candidates[idx*level+levelIdx], c, v );
				fprintf( dumpFile, "%c%d", c, v );
			}
		}
		fprintf( dumpFile, ": %d\n", support[idx]);
	}
	fprintf( dumpFile, "No. of %d node frequent episodes = %d\n\n", level, newcount);
}
