// dust.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for performing DUST filtering of query sequences

#include <stdio.h>
#include "blast.h"

// Perform dust filtering on the given ASCII sequence
void dust_dustSequence(char* originalSequence)
{
	int windowStart;
	unsigned char *sequence;
    int sequenceLength, windowLength, count;
	int cutoffScore = 20;
	int windowSize = 64;
	int minimumRegionSize = 4;
	int linker = 1;
	int windowhalf = windowSize / 2;
	struct chunk chunk;
	struct maskRegion *regions = NULL, *currentRegion = NULL;

    sequenceLength = strlen(originalSequence);

    sequence = (unsigned char*)global_malloc(sequenceLength + 1);
	strcpy(sequence, originalSequence);

    // Convert sequence into encoded format
	encoding_encodeSequence(sequence, sequenceLength, encoding_nucleotide);

	// Replace wildcards in the sequence
    count = 0;
    while (count < sequenceLength - 2)
    {
    	// If a wild
        if (sequence[count] >= encoding_numRegularLetters)
        {
            // Code replacement
            sequence[count] = encoding_randomEncodedLetter(sequence[count]);
        }
    	count++;
	}

    // Convert sequence into encoded triplets
	count = 0;
    while (count < sequenceLength - 2)
    {
        // Encode triplet
        sequence[count] = (sequence[count] << 4) | (sequence[count + 1] << 2) | sequence[count + 2];

        count++;
	}

    // Slide a window along the sequence
	for (windowStart = 0; windowStart < sequenceLength; windowStart += windowhalf)
	{
		windowLength = (int)((sequenceLength > windowStart+windowSize) ? windowSize : sequenceLength - windowStart);
		windowLength -= 2;

//        printf("process window (length=%d, position=%d)\n", windowLength, windowStart); fflush(stdout);
		dust_processWindow(windowLength, windowStart, &chunk, sequence);

//        printf("Chunk start=%d end=%d score=%d\n", chunk.start, chunk.end, chunk.score);
//    	fflush(stdout);

        // Ignore chunks that are smaller than the minimum
		if ((chunk.end - chunk.start + 1) < minimumRegionSize)
		{
			continue;
		}

		if (chunk.score > cutoffScore)
		{
        	// If this region can be linked to previous (they are close to each other)
			if (regions && regions->to + linker >= chunk.start + windowStart &&
			    regions->from <= chunk.start + windowStart)
			{
            	// Extend previous region
				regions->to = chunk.end + windowStart;
			}
			else
			{
            	// Add new region to start of linked list
				currentRegion = (struct maskRegion*)global_malloc(sizeof(struct maskRegion));
				currentRegion->from = chunk.start + windowStart;
				currentRegion->to = chunk.end + windowStart;
                currentRegion->next = regions;
                regions = currentRegion;
			}
			if (chunk.end < windowhalf)
			{
            	// Advance next window to end of chunk
				windowStart += (chunk.end - windowhalf);
			}
		}
	}

    free(sequence);

    // For each region
    currentRegion = regions;
    while (currentRegion != NULL)
    {
    	// Mask it using N's
    	count = currentRegion->from;
        while (count <= currentRegion->to)
        {
        	originalSequence[count] = 'n';
        	count++;
        }

//        printf("Start=%d End=%d\n", currentRegion->from, currentRegion->to);
        currentRegion = currentRegion->next;
    }
}

// Perform dust on the given window
void dust_processWindow(int windowLength, int windowStart, struct chunk* chunk, unsigned char* sequence)
{
	int chunkStart;

    // Initialize best chunk
	chunk->score = 0;
	chunk->start = 0;
	chunk->end = 0;

    // Get window of the sequence
    sequence += windowStart;

	// Perform dust on each chunk in the window
	for (chunkStart = 0; chunkStart < windowLength; chunkStart++)
	{
//    	printf("wo1 (%d,%d)\n", windowLength-i, i);
		dust_processChunk(windowLength - chunkStart, sequence + chunkStart, chunkStart, chunk);
	}

	// Update chunk end
	chunk->end += chunk->start;
}

// Perform dust on a chunk of a window
void dust_processChunk(int windowLength, unsigned char* sequence, int chunkStart, struct chunk* chunk)
{
    unsigned int sum;
	int position, triplet, numOccurrences;
	int newScore;
	int occurrences[256];

	// Initialize triplet occurrences to zero
    triplet = 0;
    while (triplet < 64)
    {
    	occurrences[triplet] = 0;
        triplet++;
    }

    sum = 0;
	newScore = 0;

	// For each triplet in the sequence
	for (position = 0; position < windowLength; position++)
	{
    	if (*sequence != 255)
        {
            // Increment counter of its occurance
            numOccurrences = occurrences[*sequence];

            // If it has occured more than one
            if (numOccurrences)
            {
                // Calculate score
                sum += numOccurrences;
                newScore = 10 * sum / position;

                // If the best score yet
                if (newScore > chunk->score)
                {
                    // Record the start and end of this high-scoring region
                    chunk->score = newScore;
                    chunk->start = chunkStart;
                    chunk->end = position + 2;
                }
            }
            occurrences[*sequence]++;
		}
        sequence++;
//        printf("[%d]", *occurrencesptr);
	}
}

