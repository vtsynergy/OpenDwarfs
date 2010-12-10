// readFile.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for reading the contents of a file using mmap

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

#include "blast.h"
#include "readFile.h"

// On file for reading and map to memory, then return mapped memory address
struct readFile readFile_open(char* filename)
{
	struct stat fileStats;
	struct readFile readFile;

	// Open file for reading
	if ((readFile.fileDescriptor = open(filename, O_RDONLY)) == -1)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error opening file %s for reading\n", filename);
		exit(-1);
	}

	// Get length of file
	if (fstat(readFile.fileDescriptor, &fileStats) == -1)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error opening file %s for reading\n", filename);
		exit(-1);
	}
	readFile.fileSize = fileStats.st_size;

    // Map address to fileSize bytes of application address space
    readFile.address = mmap(0, readFile.fileSize, PROT_READ, MAP_SHARED, readFile.fileDescriptor, 0);

    // Check for error in mapping
	if (readFile.address == (void*)MAP_FAILED)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error opening file %s for reading\n", filename);
		exit(-1);
	}

	return readFile;
}

// Check file exists
int readFile_checkOpen(char* filename)
{
	FILE* file;

	if ((file = fopen(filename, "r")) != NULL)
    {
    	fclose(file);
        return 1;
	}
    else
		return 0;
}

// Unmap then close the file
void readFile_close(struct readFile readFile)
{
	if (munmap(readFile.address, readFile.fileSize) < 0)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error unmapping file\n");
		exit(-1);
	}
	close(readFile.fileDescriptor);
}
