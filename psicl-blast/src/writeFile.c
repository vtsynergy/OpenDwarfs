// writeFile.c
// Copyright (c) 2005, Michael Cameron
// Permission to use this code is freely granted under the BSD license agreement,
// provided that this statement is retained.
//
// Code for writing data to file using mmap

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>

size_t writeFile_fileSize;
int writeFile_fileDescriptor;
void* writeFile_address;

// On file for reading AND writing and map to memory, then return mapped memory address
void* writeFile_open(char* filename)
{
	struct stat fileStats;

	// Open file for writing
	if ((writeFile_fileDescriptor = open(filename, O_RDWR)) == -1)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error opening file %s for writing\n", filename);
		exit(-1);
	}

	// Get length of file
	if (fstat(writeFile_fileDescriptor, &fileStats) == -1)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error getting file %s statistics\n", filename);
		exit(-1);
	}
	writeFile_fileSize = fileStats.st_size;

	// Map address to fileSize bytes of application address space
	writeFile_address = mmap(0, writeFile_fileSize, PROT_WRITE, MAP_SHARED, writeFile_fileDescriptor, 0);

	// Check for error in mapping
	if (writeFile_address == (void*)MAP_FAILED)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error mapping file %s data to program address space\n", filename);
		exit(-1);
	}

	return writeFile_address;
}

// Get file size/length
size_t writeFile_getSize()
{
	return writeFile_fileSize;
}

// Unmap then close the file
void writeFile_close()
{
	if (munmap(writeFile_address, writeFile_fileSize) < 0)
	{
        fprintf(stderr, "%s\n", strerror(errno));
		fprintf(stderr, "Error unmapping file\n");
		exit(-1);
	}
	close(writeFile_fileDescriptor);
}
