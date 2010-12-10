// createindex.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Creates an .index file for a database for inverted index mode search

#include "blast.h"
#include <errno.h>
#include <stdio.h>

int4 main(int4 argc, char* argv[])
{
	unsigned char *filename, *address, *sequence, code, *startAddress, remainder;
	uint4 numberOfSequences, longestSequenceLength, dbAlphabetType;
    uint8 numberOfLetters;
	uint4 descriptionStart = 0, descriptionLength = 0, sequenceSize, count;
	uint4 fileSize, sequenceCount = 0, encodedLength;
	uint4 allocSize = 1024;
	uint4* int4Address;
	char* sequenceFilename, *indexFilename;
	struct readFile readFile;
    uint4 *sequencePositions, *descriptionLocations;
    uint4 descriptionFileSize = 0;

    uint4 numWords, codeword, numOffsets, indexSize = 0;
	unsigned char* offsets, *vbyteEncoded, startVbyteEncoded[4];
	FILE *indexFile;
	unsigned char headerData[40], *headerDataPointer;
    uint4 headerSize, positionBlockOffset = 0;
	uint4* listPositions, currentPosition = 0;
	uint4 fromCodeword, toCodeword, codewordRangeSize;

    // User must provide FASTA format file at command line
	if (argc != 2 && argc != 4)
	{
		fprintf(stderr, "Useage: createindex <Database> [Word size] [Interval size]\n");
		exit(-1);
	}
	filename = argv[1];

    // Optionally user provides word and interval size
    if (argc == 4)
    {
    	index_wordSize = atoi(argv[2]);
    	index_intervalSize = atoi(argv[3]);
    }

	// Allocation initial buffer to hold sequences
	sequence = (char*)global_malloc(allocSize);

	// Construct sequence filename
	sequenceFilename = (char*)global_malloc(strlen(filename) + 13);
	sprintf(sequenceFilename, "%s.sequences", filename);

	// Open sequence file for reading, mapping contents to address
	readFile = readFile_open(sequenceFilename);
	startAddress = address = (char*)readFile.address;

    // Get the file size
	fileSize = readFile.fileSize;

	// Read database statistics (First 40 bytes contains 4 vbytes)
    vbyte_getVbyte(address, &numberOfSequences);
    address = vbyte_get64vbyte(address, &numberOfLetters);
    vbyte_getVbyte(address, &longestSequenceLength);
    vbyte_getVbyte(address, &dbAlphabetType);

    // Advance to start of sequences
    address = startAddress + 40;

    // Initialize list of sequence and description positions
    sequencePositions = (uint4*)malloc(sizeof(uint4) * numberOfSequences);
    descriptionLocations = (uint4*)malloc(sizeof(uint4) * numberOfSequences);

    // Initialize codes array
	encoding_initialize(dbAlphabetType);

    printf("Processing database...");
	fflush(stdout);

    // Move through sequence file reading sequence codes and converting back to ASCII sequences
	while (address < startAddress + fileSize)
	{
        // Record start position of sequence
		sequencePositions[sequenceCount] = address - startAddress;

        // Read two vbytes, the description location in the FASTA file and the length of the sequence
        descriptionStart += descriptionLength;
        vbyte_getVbyte(address, &descriptionLength);
        vbyte_getVbyte(address, &sequenceSize);

        // Record location of description
        descriptionLocations[sequenceCount] = descriptionStart;

        // Nucleotide sequences
        if (encoding_alphabetType == encoding_nucleotide)
        {
            // Read a third vbyte; the total length of sequence data including wildcard info
            vbyte_getVbyte(address, &encodedLength);

            // Skip to end of the sequence
            address += encodedLength;
        }
        // Protein sequences
        else
        {
            // Skip past sentinal byte
            address++;

            sequence = address;
            address += sequenceSize;

            // Skip past sentinal byte
            address++;
		}

		sequenceCount++;

        // Print status dots
		if (sequenceCount % 10000 == 0)
		{
			printf(".");
			fflush(stdout);
		}
	}

	// Check we ended at the end of the file
	if (address != startAddress + fileSize)
	{
		fprintf(stderr, "Error: Premature end of sequence file\n");
		exit(-1);
	}

	printf("done.\n");
	printf("%d sequences read.\n", sequenceCount);
	fflush(stdout);

    // Create the index
	printf("Creating index...");
	fflush(stdout);

    indexFilename = (char*)global_malloc(strlen(filename) + 10);
	sprintf(indexFilename, "%s.index", filename);

	// Open index file for writing
	if ((indexFile = fopen(indexFilename, "w")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for writing\n", indexFilename);
		exit(-1);
	}

    // Write word size and interval size at start of index
    headerDataPointer = headerData;
    vbyte_putVbyte(headerDataPointer, index_wordSize);
    vbyte_putVbyte(headerDataPointer, index_intervalSize);
    headerSize = headerDataPointer - headerData;
    if (fwrite(&headerData, sizeof(unsigned char), headerSize, indexFile) < headerSize)
	{
		fprintf(stderr, "Error writing header to index file %s\n", indexFilename);
		exit(-1);
	}

    // Write sequence positions
    if (fwrite(sequencePositions, sizeof(uint4), numberOfSequences, indexFile) < numberOfSequences)
    {
        fprintf(stderr, "Error writing sequence positions to index file %s\n", indexFilename);
    	fprintf(stderr, strerror(errno));
        fprintf(stderr, "\n"); fflush(stderr);
        exit(-1);
    }

    // Write description locations
    if (fwrite(descriptionLocations, sizeof(uint4), numberOfSequences, indexFile) < numberOfSequences)
    {
        fprintf(stderr, "Error writing description locations to index file %s\n", indexFilename);
    	fprintf(stderr, strerror(errno));
        fprintf(stderr, "\n"); fflush(stderr);
        exit(-1);
    }

    // Write padding for list positions
    numWords = pow(4, index_wordSize);
	listPositions = (uint4*)global_malloc(sizeof(uint4) * (numWords + 1));
	codeword = 0;
    while (codeword < numWords + 1)
    {
		listPositions[codeword] = 0;
    	codeword++;
	}

    if (fwrite(listPositions, sizeof(uint4), numWords + 1, indexFile) < numWords + 1)
    {
        fprintf(stderr, "Error writing word offset position padding to index file %s\n", indexFilename);
    	fprintf(stderr, strerror(errno));
        fprintf(stderr, "\n"); fflush(stderr);
        exit(-1);
    }

    // SECOND PASS: build inverted index lists

    // Determine the number of words to process each pass through the collection
    codewordRangeSize = (300000000 / (sizeof(struct wordList) +
                        (6.0 * (float)numberOfLetters / (float)numWords / (float)index_intervalSize)));
	printf("codewordRangeSize=%d\n", codewordRangeSize);

	// Pass through the collection considering words in the range between fromCodeword and toCodeword
    fromCodeword = 0;
    toCodeword = 0;
    while (toCodeword < numWords)
    {
    	// Set fromCodeword and toCodeword
		fromCodeword = toCodeword;
        toCodeword += codewordRangeSize;
        if (toCodeword > numWords)
        	toCodeword = numWords;

		printf("[%d,%d]/%d\n", fromCodeword, toCodeword, numWords);

        // Initialize index build
        index_initializeBuild(fromCodeword, toCodeword);

        // Advance to start of sequences
        address = startAddress + 40;

        // Move through sequence file reading sequences
        while (address < startAddress + fileSize)
        {
            // Read two vbytes, the description location in the FASTA file and the length of the sequence
            descriptionStart += descriptionLength;
            vbyte_getVbyte(address, &descriptionLength);
            vbyte_getVbyte(address, &sequenceSize);

            // Nucleotide sequences
            if (encoding_alphabetType == encoding_nucleotide)
            {
                // Read a third vbyte; the total length of sequence data including wildcard info
                vbyte_getVbyte(address, &encodedLength);

                // Unpack the sequence
                sequence = encoding_byteUnpack(address, sequenceSize);

                // Skip to end of the sequence
                address += encodedLength;
            }
            // Protein sequences
            else
            {
                // Skip past sentinal byte
                address++;

                sequence = address;
                address += sequenceSize;

                // Skip past sentinal byte
                address++;
            }

            // Add the subject to the index structure
            index_addSubject(sequence, sequenceSize, fromCodeword, toCodeword);

            if (encoding_alphabetType == encoding_nucleotide)
                free(sequence);

            sequenceCount++;

            // Print status dots
            if (sequenceCount % 10000 == 0)
            {
                printf(".");
                fflush(stdout);
            }
        }

        // Check we ended at the end of the file
        if (address != startAddress + fileSize)
        {
            fprintf(stderr, "Error: Premature end of sequence file\n");
            exit(-1);
        }

        // Get word list positions and number of words
//        wordOffsetPositions = index_wordOffsetPositions(fromCodeword, toCodeword);
//        numWords = toCodeword - fromCodeword;

        // For each word, write list of offset-gaps to END of file
        codeword = fromCodeword;
        while (codeword < toCodeword)
        {
            // Write the list of offsets
            numOffsets = index_numWordOffsets(codeword);
            offsets = index_wordOffsets(codeword);

            // Record position of offsets
			listPositions[codeword] = currentPosition;
            currentPosition += numOffsets;

            if (numOffsets > 0)
            {
/*                int a = 0, b = 0;
                while (a < numOffsets)
                {
                    b += offsets[a];
                    a++;
                }
                fprintf(stderr, "[%p][%d][%p][%d]\n", offsets, numOffsets, indexFile, b); fflush(stderr);
*/
//                if (fwrite("abc", sizeof(unsigned char), 3, indexFile) < 3)
                if (fwrite(offsets, sizeof(unsigned char), numOffsets, indexFile) < numOffsets)
                {
                    fprintf(stderr, "Error writing word offsets to index file %s\n", indexFilename);
                    fprintf(stderr, strerror(errno));
                    fprintf(stderr, "\n"); fflush(stderr);
                    exit(-1);
                }
			}

            indexSize += numOffsets;
            codeword++;
        }

        // Finishing building inverted lists for the given range of words
        index_finishBuild(fromCodeword, toCodeword);
    }

    readFile_close(readFile);

    // Record extra list position
    listPositions[toCodeword] = currentPosition;

    // Close writing to index file
    fclose(indexFile);

    // Reopen and write list positions
	if ((indexFile = fopen(indexFilename, "r+")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for writing\n", indexFilename);
		exit(-1);
	}

    // Jump to position of interest
    fseek(indexFile, numberOfSequences * sizeof(uint4) * 2 + headerSize, SEEK_SET);

    // Write list positions to index file
    if (fwrite(listPositions, sizeof(uint4), numWords + 1, indexFile) < numWords + 1)
    {
        fprintf(stderr, "Error writing word offset positions to index file %s\n", indexFilename);
    	fprintf(stderr, strerror(errno));
        fprintf(stderr, "\n"); fflush(stderr);
        exit(-1);
    }

	// Close writing to index file
	fclose(indexFile);

    printf("done.\n");

    printf("Size of index = %dMb + %dMb + %dMb + %dMb = ", numberOfSequences * 4 / 1048576,
           numberOfSequences * 4 / 1048576, ((numWords + 1) * 4) / 1048576,
           indexSize / 1048576);

    indexSize += (numWords + 1) * 4 + numberOfSequences * 4 + numberOfSequences * 4;
    printf("%dMb\n", indexSize / 1048576);

	return 0;
}
