// dustApp
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Standalone tool for DUST filtering sequences

#include <stdio.h>
#include "blast.h"

int main(int argc, char* argv[])
{
	char *filename, *sequence, *description, *sequenceCopy;
    int sequenceLength;

	// User must provide FASTA format file at command line
	if (argc < 2)
	{
		fprintf(stderr, "Useage: dust <FASTA file>\n");
		exit(-1);
	}
	filename = argv[1];

    // Initialize encoding routines
    encoding_initialize(encoding_nucleotide);

    // Open FASTA file for reading
	readFasta_open(filename);

    // Read each sequence from the file
    while (readFasta_readSequence())
	{
		// Get sequence just read
		sequence = readFasta_sequenceBuffer;
		description = readFasta_descriptionBuffer;
        sequenceLength = readFasta_sequenceLength;

        // Make in-memory copy of it
        sequenceCopy = (char*)global_malloc(sequenceLength);
        strcpy(sequenceCopy, sequence);

		// Perform dust filtering
        dust_dustSequence(sequence);

        // Print description and filtering sequence
        printf(">%s\n%s\n", description, sequence);
    }
}

