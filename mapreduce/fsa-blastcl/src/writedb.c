// writedb.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Functions for creating formatted BLAST .sequences and .descriptions files

#include "blast.h"

char* writedb_filename;
FILE *writedb_sequenceFile, *writedb_descriptionsFile, *writedb_dataFile;
uint8 writedb_volumeSize;
char *writedb_sequenceFilename, *writedb_descriptionsFilename, *writedb_dataFilename;
uint4 writedb_maximumSequenceLength, writedb_alphabetType, writedb_minimumSequenceLength;
uint8 writedb_numberOfLetters;
uint4 writedb_volume, writedb_sequenceCount, writedb_numberOfClusters;
struct memBlocks* writedb_sequenceData;

// Initialize writing to formatted database
void writedb_initialize(char* filename, uint4 alphabetType)
{
	char* wildcardsFilename;

    writedb_filename = filename;
    writedb_alphabetType = alphabetType;
	writedb_maximumSequenceLength = 0;
	writedb_minimumSequenceLength = 0;
	writedb_numberOfLetters = 0;
	writedb_volume = 0;
    writedb_sequenceCount = 0;
    writedb_numberOfClusters = 0;

	// Construct sequence and description filenames
	writedb_sequenceFilename = (char*)global_malloc(strlen(filename) + 13);
	sprintf(writedb_sequenceFilename, "%s.sequences", filename);
	writedb_descriptionsFilename = (char*)global_malloc(strlen(filename) + 15);
	sprintf(writedb_descriptionsFilename, "%s.descriptions", filename);
	writedb_dataFilename = (char*)global_malloc(strlen(filename) + 8);
	sprintf(writedb_dataFilename, "%s.data", filename);
	wildcardsFilename = (char*)global_malloc(strlen(filename) + 12);
	sprintf(wildcardsFilename, "%s.wildcards", filename);

    // Delete the wildcards file if one exists
    rename(wildcardsFilename, writedb_sequenceFilename);

	// Open sequence file for writing
	if ((writedb_sequenceFile = fopen(writedb_sequenceFilename, "w")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for writing\n", writedb_sequenceFilename);
		exit(-1);
	}

    // Write sentinal/padding byte at start
	if (alphabetType == encoding_protein)
        fputc(encoding_sentinalCode, writedb_sequenceFile);
    else
        fputc(0, writedb_sequenceFile);

	// Open descriptions file for writing
	if ((writedb_descriptionsFile = fopen(writedb_descriptionsFilename, "w")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for writing\n", writedb_descriptionsFilename);
		exit(-1);
	}

    writedb_volumeSize = 1;
    writedb_sequenceData = memBlocks_initialize(sizeof(struct sequenceData),
                           constants_initialSequenceData);
}

// Add sequence to the formatted collection
void writedb_addSequence(unsigned char* sequence, uint4 sequenceLength, unsigned char* description,
                         uint4 descriptionLength, unsigned char* wildcards, uint4 wildcardsLength,
                         struct child* children, uint4 numChildren)
{
	uint4 encodedLength, childNum, sizeEdits = 0, editNum;
	unsigned char *editData, *startEditData;
    struct child child;
    struct sequenceData* sequenceData;

    sequenceData = memBlocks_newEntry(writedb_sequenceData);

    // Write the description to file
    if (description != NULL)
        if (fwrite(description, sizeof(unsigned char), descriptionLength, writedb_descriptionsFile)
            < descriptionLength)
        {
            fprintf(stderr, "Error writing header to sequence file %s\n", writedb_sequenceFilename);
            exit(-1);
        }

    // Calculate length of encoded sequence
    if (writedb_alphabetType == encoding_nucleotide)
    {
        encodedLength = encoding_bytePackSequence(sequence, sequenceLength);
	}
    else
    {
    	encodedLength = sequenceLength + 2;
    }

	// Calculate maximum space required to record sequence's edits
    childNum = 0;
	while (childNum < numChildren)
    {
    	child = children[childNum];
		sizeEdits += 16 + 5 * child.numEdits;
        childNum++;
	}

    // Initialize array to record edits
   	editData = startEditData = global_malloc(sizeEdits);

    // Record children edits as vbytes
    childNum = 0;
	while (childNum < numChildren)
    {
    	child = children[childNum];

        // Write children descriptions to disk
        if (fwrite(child.description, sizeof(unsigned char), child.descriptionLength,
            writedb_descriptionsFile) < child.descriptionLength)
        {
            fprintf(stderr, "Error writing description to sequence file %s\n", writedb_descriptionsFilename);
            exit(-1);
        }
        descriptionLength += child.descriptionLength;

        // Convert child details to vbytes
        vbyte_safePutVbyte(editData, child.descriptionLength);
        vbyte_safePutVbyte(editData, child.regionStart);
        vbyte_safePutVbyte(editData, child.length);
        vbyte_safePutVbyte(editData, child.numEdits);

        // Append edits
        editNum = 0;
        while (editNum < child.numEdits)
        {
        	// Record edit character
			*editData = child.edits[editNum].code;
            editData++;

        	editNum++;
        }

        // Add sequence size to total tally of letters
        writedb_numberOfLetters += child.length;
        writedb_sequenceCount++;

        childNum++;
    }

    // Update volume size, encoded length
    encodedLength += (editData - startEditData);
    writedb_volumeSize += encodedLength + wildcardsLength;

    sequenceData->descriptionLength = descriptionLength;
	sequenceData->sequenceLength = sequenceLength;
	sequenceData->encodedLength = encodedLength + wildcardsLength;

    // If the entry will exceed volume max size
    if (writedb_volumeSize > constants_volumeMaxSize)
    {
        // Close current volume
        fclose(writedb_sequenceFile);

        // Open next volume for writing
        writedb_volume++;
        sprintf(writedb_sequenceFilename, "%s.sequences%d", writedb_filename, writedb_volume);
        if ((writedb_sequenceFile = fopen(writedb_sequenceFilename, "w")) == NULL)
        {
            fprintf(stderr, "Error opening file %s for writing\n", writedb_sequenceFilename);
            exit(-1);
        }

        // Reset volume size counter
        writedb_volumeSize = encodedLength + wildcardsLength;
    }

    // Nulceotide
    if (writedb_alphabetType == encoding_nucleotide)
    {
        // Write packed nucleotide sequences to disk
        if (fwrite(sequence, sizeof(unsigned char), encodedLength, writedb_sequenceFile) < encodedLength)
        {
            fprintf(stderr, "Error writing to sequence file %s\n", writedb_sequenceFilename);
            exit(-1);
        }
    }
    // Protein
    else
    {
        // Write sentinal byte after protein sequences
        fputc(encoding_sentinalCode, writedb_sequenceFile);

        // Write sequence codes to disk
        if (fwrite(sequence, sizeof(unsigned char), sequenceLength, writedb_sequenceFile) < sequenceLength)
        {
            fprintf(stderr, "Error writing to sequence file %s\n", writedb_sequenceFilename);
            exit(-1);
        }

        // Write sentinal byte after protein sequences
        fputc(encoding_sentinalCode, writedb_sequenceFile);
    }

    // Write wildcard data to disk
    if (fwrite(wildcards, sizeof(unsigned char), wildcardsLength, writedb_sequenceFile) < wildcardsLength)
    {
        fprintf(stderr, "Error writing to sequence file %s\n", writedb_sequenceFilename);
        exit(-1);
    }

    // Write edit information to disk
    if (fwrite(startEditData, sizeof(unsigned char), (editData - startEditData),
               writedb_sequenceFile) < (editData - startEditData))
    {
        fprintf(stderr, "Error writing to sequence file %s\n", writedb_sequenceFilename);
        exit(-1);
    }
    free(startEditData);

	if (numChildren == 0)
    {
        // Add sequence size to total tally of letters
        writedb_numberOfLetters += sequenceLength;
        writedb_sequenceCount++;
	}

    writedb_numberOfClusters++;

    // Check for new longest/shortest sequence
    if (sequenceLength > writedb_maximumSequenceLength)
        writedb_maximumSequenceLength = sequenceLength;
    if (writedb_minimumSequenceLength == 0 || sequenceLength < writedb_minimumSequenceLength)
        writedb_minimumSequenceLength = sequenceLength;
}

// Finalize writing to the formatted collection
void writedb_close()
{
	unsigned char headerData[40], *headerDataPointer;
	uint4 headerLength;
    struct sequenceData* sequenceData;

    // Write sentinal/padding byte at end
	if (writedb_alphabetType == encoding_protein)
        fputc(encoding_sentinalCode, writedb_sequenceFile);
    else
        fputc(0, writedb_sequenceFile);

    // Close writing to sequence and description files
	fclose(writedb_sequenceFile);
	fclose(writedb_descriptionsFile);

	// Open data file for writing
	if ((writedb_dataFile = fopen(writedb_dataFilename, "w")) == NULL)
	{
		fprintf(stderr, "Error opening file %s for writing\n", writedb_dataFilename);
		exit(-1);
	}

	// Convert 6 header values to vbytes
    headerDataPointer = headerData;
    vbyte_safePutVbyte(headerDataPointer, constants_databaseVersion);
    vbyte_safePutVbyte(headerDataPointer, writedb_sequenceCount);
    vbyte_safePutVbyte(headerDataPointer, writedb_numberOfLetters);
    vbyte_safePutVbyte(headerDataPointer, writedb_maximumSequenceLength);
    vbyte_safePutVbyte(headerDataPointer, writedb_alphabetType);
    vbyte_safePutVbyte(headerDataPointer, writedb_numberOfClusters);
	vbyte_safePutVbyte(headerDataPointer, writedb_volume + 1);

    // Write the header data at the start of the file
    headerLength = headerDataPointer - headerData;
    if (fwrite(&headerData, sizeof(unsigned char), headerLength, writedb_dataFile) < headerLength)
	{
		fprintf(stderr, "Error writing header to sequence file %s\n", writedb_dataFilename);
		exit(-1);
	}

    // For each sequence
	memBlocks_resetCurrent(writedb_sequenceData);
    while ((sequenceData = memBlocks_getCurrent(writedb_sequenceData)) != NULL)
    {
        // Prepare to write sequence description length, subject length, and encoded length using vbytes
	    headerDataPointer = headerData;
        vbyte_safePutVbyte(headerDataPointer, sequenceData->descriptionLength);
        vbyte_safePutVbyte(headerDataPointer, sequenceData->sequenceLength);
        vbyte_safePutVbyte(headerDataPointer, sequenceData->encodedLength);

        // Write sequence header information
        headerLength = headerDataPointer - headerData;
        if (fwrite(headerData, sizeof(unsigned char), headerLength, writedb_dataFile) < headerLength)
        {
            fprintf(stderr, "Error writing to sequence file %s\n", writedb_dataFilename);
            exit(-1);
        }
	}

	// Close writing to sequence file
	fclose(writedb_dataFile);

    memBlocks_free(writedb_sequenceData);
}

