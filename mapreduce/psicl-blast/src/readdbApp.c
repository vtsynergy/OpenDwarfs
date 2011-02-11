// readdbApp.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Debugging program that reads .sequences and .descriptions files and prints their contents

#include "blast.h"

int4 main(int4 argc, char* argv[])
{
	unsigned char *filename, *sequence;
	uint4 descriptionStart = 0, descriptionLength = 0, sequenceLength;
	uint4 encodedLength, numChildren, childNum, count;
	char *description;
    struct child* children, *child;
	uint4* clusterSizes = NULL, numClusterSizes = 0;
	uint4 display = 1;

	// User must provide FASTA format file at command line
	if (argc < 2)
	{
		fprintf(stderr, "Useage: readdb <FASTA file>\n");
		exit(-1);
	}
	filename = argv[1];

	readdb_open(filename);

    printf("Number of clusters = %u\n", readdb_numberOfClusters);
    printf("Number of sequences = %u\n", readdb_numberOfSequences);
    printf("Number of volumes = %u\n", readdb_numberOfVolumes);
	printf("Total number of letters = %llu\n", readdb_numberOfLetters);
	printf("Length of longest sequence = %u\n", readdb_longestSequenceLength);
	printf("Alphabet type = %s\n", encoding_alphabetTypes[readdb_dbAlphabetType]);

	// Initialize codes array
	encoding_initialize(readdb_dbAlphabetType);

    do
    {
        // Read each sequence in the collection
        while (readdb_readSequence(&sequence, &sequenceLength, &descriptionStart,
                                   &descriptionLength, &encodedLength))
        {
            // Unpack nucleotide sequences
            if (encoding_alphabetType == encoding_nucleotide)
                sequence = encoding_byteUnpack(sequence, sequenceLength);

            if (encoding_alphabetType == encoding_protein && sequenceLength + 2 != encodedLength)
            {
                // Get the children
                children = readdb_getChildren(sequence, sequenceLength, encodedLength,
                                              descriptionStart, &numChildren);

				// Record number of clusters of each size
				if (numChildren + 1 > numClusterSizes)
                {
					clusterSizes = global_realloc(clusterSizes, sizeof(uint4) * (numChildren + 1));
                    while (numClusterSizes < numChildren + 1)
                    {
						clusterSizes[numClusterSizes] = 0;
                    	numClusterSizes++;
                    }
                }
                clusterSizes[numChildren]++;

//                if (sequenceLength < 60)
//                display = 1;

                if (display)
                printf("\n*** Parent with %d children ***\n", numChildren);

                // Print cluster
                if (display)
                {
                	print_singleSequence(sequence, sequenceLength); printf("\n");
                }

                // For each child
                childNum = 0;
                while (childNum < numChildren)
                {
                    child = children + childNum;

                    // Align with parent
                    count = 0;
                    if (display)
                    while (count < child->regionStart)
                    {
                        printf(" ");
                        count++;
                    }

                    // Print child sequence and description
                    if (display)
                    print_singleSequence(child->sequence, child->length);

                    description = descriptions_getDescription(child->descriptionLocation,
                                                              child->descriptionLength);

                    // Print descriptions at end
                    count = child->regionStart + child->length;
                    if (display)
                    while (count < sequenceLength)
                    {
                        printf(" ");
                        count++;
                    }

                    if (display)
                    {
//                    	printf(" (%s dloc=%d)\n", description, child->descriptionLocation);
                    	printf(" (%s)\n", description);
					}

                    free(child->sequence - 1);
                    free(child->edits);

                    childNum++;
                }

                free(children);
            }
            else
            {
                // Print the sequence description
                description = descriptions_getDescription(descriptionStart, descriptionLength);
                if (display)
                {
//                	printf(">%s (dloc = %d)\n", description, descriptionStart);
                	printf(">%s\n", description);
				}

                // Print sequence
                if (display)
                {
                	print_singleSequence(sequence, sequenceLength); printf("\n");
                }
            }

            // Free unpacked sequence
            if (encoding_alphabetType == encoding_nucleotide)
                free(sequence);
        }
	}
    while (readdb_nextVolume());

	printf("%d sequences read.\n", readdb_numberOfSequences);
	fflush(stdout);

    while (numClusterSizes > 0)
    {
    	numClusterSizes--;
        printf("%d clusters with %d children\n", clusterSizes[numClusterSizes], numClusterSizes);
    }

	return 0;
}
