// formatdb.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Program for creating a database .sequences file from a FASTA format database

#include <blast.h>

uint4 determineDbAlphabetType(char* filename);

int4 main(int argc, char* argv[])
{
	char *sequence, *filename;
    uint4 sequenceLength;
    int4 totalWilds = 0, alphabetType;
    struct memSingleBlock* wildcardEdits;
    struct wildcardEdit* wildcardEdit;
    char *wildcardData = NULL, *startWildcardData = NULL;

	// User must provide FASTA format file at command line
	if (argc < 2)
	{
		fprintf(stderr, "Useage: formatdb <FASTA file>\n");
		exit(-1);
	}
	filename = argv[1];

    // Initialize array to store wildcard edits
    wildcardEdits = memSingleBlock_initialize(sizeof(struct wildcardEdit), 10);

    // Determine if database is protein or nucleotide
	alphabetType = determineDbAlphabetType(filename);

    if (alphabetType == encoding_protein)
    {
    	printf("PROTEIN database detected.\n");
    }
    else if (alphabetType == encoding_nucleotide)
    {
    	printf("NUCLEOTIDE database detected.\n");
    }

	// Initialize codes array
	encoding_initialize(alphabetType);

	// Initialize writing to formatted database
    writedb_initialize(filename, alphabetType);

    // Open FASTA file for reading
	readFasta_open(filename);

	printf("Formatting database...");
	fflush(stdout);

	// Move through the FASTA file reading descriptions and sequences
	while (readFasta_readSequence())
	{
		// Get sequence just read
		sequence = readFasta_sequenceBuffer;
		sequenceLength = readFasta_sequenceLength;

        // Encode the sequence
		encoding_encodeSequence(sequence, sequenceLength, alphabetType);

        // Convert nucleotide sequences to byte-packed format
        if (alphabetType == encoding_nucleotide)
        {
			// Replace any wilds with a random character
            totalWilds += encoding_replaceWildcards(wildcardEdits, sequence, sequenceLength);

			// Declare memory to hold wildcard data
            startWildcardData = global_realloc(startWildcardData,
                                sizeof(char) * wildcardEdits->numEntries * 5);
            wildcardData = startWildcardData;

            // For each wildcard edit, encode details using chars and vbytes
            memSingleBlock_resetCurrent(wildcardEdits);
            while ((wildcardEdit = memSingleBlock_getCurrent(wildcardEdits)) != NULL)
            {
                // Record wild character
                *wildcardData = wildcardEdit->code;
                wildcardData++;

                // Convert the position to a vbyte
                vbyte_putVbyte(wildcardData, wildcardEdit->position);
            }
		}
        else
        {
        	startWildcardData = wildcardData = NULL;
		}

//        printf("[%s](%d)", readFasta_descriptionBuffer, readFasta_descriptionLength); fflush(stdout);

        // Add sequence to the formatted collection
        writedb_addSequence(sequence, sequenceLength, readFasta_descriptionBuffer,
                            readFasta_descriptionLength, startWildcardData,
                            wildcardData - startWildcardData, NULL, 0);

		// Print status dots
		if (writedb_sequenceCount % 10000 == 0)
		{
			printf(".");
			fflush(stdout);
		}
	}

	// Close fasta reader
	readFasta_close();

    // Finalize writing to the formatted collection
    writedb_close();

	printf("done.\n");
	printf("%d sequences processed.\n", writedb_sequenceCount);
	printf("%llu letters processed.\n", writedb_numberOfLetters);
    printf("%d wildcards encoded.\n", totalWilds);
	printf("%d volume(s) created.\n", writedb_volume + 1);
	printf("Longest/shortest sequence was %d/%d letters\n",
           writedb_maximumSequenceLength, writedb_minimumSequenceLength);
	fflush(stdout);

	return 0;
}

// Read the first 10 sequences from the database to determine its type
uint4 determineDbAlphabetType(char* filename)
{
	int4 sequenceCount = 0;
    char* sequence;
    uint4 sequenceLength;

    // Open FASTA file for reading
	readFasta_open(filename);

	// Move through the FASTA file reading descriptions and sequences
	while (readFasta_readSequence() && sequenceCount < 10)
	{
		// Get sequence just read
		sequence = readFasta_sequenceBuffer;
		sequenceLength = readFasta_sequenceLength;

        // Determine the alphabet of the current sequence
        if (encoding_determineAlphabetType(sequence, sequenceLength) == encoding_protein)
		{
        	// If contains protein letters, return protein type
            readFasta_close();
			return encoding_protein;
        }

        sequenceCount++;
	}

	// Close fasta reader and return nucleotide type
	readFasta_close();
    return encoding_nucleotide;
}
