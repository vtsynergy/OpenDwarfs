// readFasta.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for reading a FASTA format file and extracting sequence and description
// information

#include "blast.h"
#include <stdio.h>
#include <errno.h>

char* readFasta_descriptionBuffer;
char* readFasta_sequenceBuffer;
uint4 readFasta_sequenceLength;
uint4 readFasta_descriptionLength;

int4 readFasta_bufferAlloc;
int4 readFasta_sequenceBufferAlloc;
char* readFasta_lineBuffer;
char readFasta_endReached;

char* readFasta_filename;
FILE* readFasta_file;
int4 readFasta_lineLength;

// Open FASTA file for reading
void readFasta_open(char* filename)
{
	readFasta_file = fopen(filename, "r");

    if (readFasta_file == NULL)
    {
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error opening file %s for reading\n", filename);
		exit(-1);
    }

    readFasta_filename = filename;

	readFasta_bufferAlloc = 10;
	readFasta_sequenceBufferAlloc = 10;

    // Declare memory for buffers
	readFasta_sequenceBuffer = (char*)global_malloc(sizeof(char) * readFasta_bufferAlloc);
	readFasta_descriptionBuffer = (char*)global_malloc(sizeof(char) * readFasta_bufferAlloc);
	readFasta_lineBuffer = (char*)global_malloc(sizeof(char) * readFasta_bufferAlloc);
    readFasta_lineBuffer[0] = '\0'; readFasta_lineBuffer[1] = '\0';
    readFasta_lineLength = 1;
    readFasta_endReached = 0;
}

// Read a description and sequence from file, return 0 if end of file, otherwise 1
int4 readFasta_readSequence()
{
	uint4 pos, pos2;

	// We've reached end of file
	if (readFasta_endReached)
    	return 0;

	// Clear sequence buffer
	readFasta_sequenceBuffer[0] = '\0';
    readFasta_sequenceLength = 0;

    // Line buffer contains description for new sequence
	memcpy(readFasta_descriptionBuffer, readFasta_lineBuffer + 1, readFasta_lineLength);
    readFasta_descriptionLength = readFasta_lineLength - 1;

	// Read each line from FASTA file
	while (!readFasta_endReached && fgets(readFasta_lineBuffer, readFasta_bufferAlloc, readFasta_file))
    {
    	// If we didn't read an entire line
        readFasta_lineLength = strlen(readFasta_lineBuffer);
    	while (readFasta_lineBuffer[readFasta_lineLength - 1] != '\n')
        {
//        	printf("Read bit \"%s\"\n", readFasta_lineBuffer);

            // Double size of the buffers and read more
			readFasta_lineBuffer = (char*)global_realloc(readFasta_lineBuffer,
            	                   sizeof(char) * readFasta_bufferAlloc * 2);
			readFasta_descriptionBuffer = (char*)global_realloc(readFasta_descriptionBuffer,
            	                          sizeof(char) * readFasta_bufferAlloc * 2);

			// Read next line
			if (!fgets(readFasta_lineBuffer + readFasta_bufferAlloc - 1, readFasta_bufferAlloc, readFasta_file))
			{
            	// If reached end of file
                readFasta_lineLength = strlen(readFasta_lineBuffer);
                readFasta_endReached = 1;
                break;
			}

            readFasta_bufferAlloc = readFasta_bufferAlloc * 2 - 1;
            readFasta_lineLength = strlen(readFasta_lineBuffer);
        }

        // Remove the trailing \n
        if (readFasta_lineBuffer[readFasta_lineLength - 1] == '\n')
        {
            readFasta_lineLength--;
            readFasta_lineBuffer[readFasta_lineLength] = '\0';
		}

//        printf("Read line \"%s\"\n", readFasta_lineBuffer);

        // If this is a description line
        if (readFasta_lineBuffer[0] == '>')
        {
//        	printf("DESCRIPTION!\n");
			// If the description buffer is not empty
            if (readFasta_descriptionBuffer[0] != '\0')
			{
            	return 1;
            }

//        	printf("(First)\n");
            // Else line is description of first sequence
            memcpy(readFasta_descriptionBuffer, readFasta_lineBuffer + 1, readFasta_lineLength);
            readFasta_descriptionLength = readFasta_lineLength - 1;

//            printf("(end copy)\n");
        }
        // Otherwise it is part of a sequence
        else
        {
        	// If there is insufficient space in the current sequence buffer
			if (readFasta_lineLength + readFasta_sequenceLength + 1 > readFasta_sequenceBufferAlloc)
            {
                // Increase size of the buffer
				readFasta_sequenceBufferAlloc = (readFasta_lineLength + readFasta_sequenceLength + 1) * 2;
                readFasta_sequenceBuffer = (char*)global_realloc(readFasta_sequenceBuffer,
                                           sizeof(char) * readFasta_sequenceBufferAlloc);
			}

            // Remove non-alphabet characters from the sequence
			pos = 0; pos2 = 0;
            while (pos < readFasta_lineLength)
            {
            	if ((readFasta_lineBuffer[pos] >= 'a' && readFasta_lineBuffer[pos] <= 'z') ||
                    (readFasta_lineBuffer[pos] >= 'A' && readFasta_lineBuffer[pos] <= 'Z'))
                {
                	readFasta_lineBuffer[pos2] = readFasta_lineBuffer[pos];
                    pos2++;
                }
            	pos++;
            }
			readFasta_lineBuffer[pos2] = '\0';
			readFasta_lineLength = pos2;

            // Append new part of sequence onto the current sequence buffer contents
			memcpy(readFasta_sequenceBuffer + readFasta_sequenceLength, readFasta_lineBuffer,
                   readFasta_lineLength + 1);
            readFasta_sequenceLength += readFasta_lineLength;
        }
    }

    // End of file. Return contents of last sequence
    readFasta_endReached = 1;
    return 1;
}

// Close the FASTA file and free memory
void readFasta_close()
{
	free(readFasta_lineBuffer);
	free(readFasta_sequenceBuffer);
	free(readFasta_descriptionBuffer);
    fclose(readFasta_file);
}
