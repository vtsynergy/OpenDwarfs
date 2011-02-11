#ifndef _writeFile_
#define _writeFile_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

// On file for reading AND writing and map to memory, then return mapped memory address
void* writeFile_open(char* filename)

// Get file size/length
int4 writeFile_getSize()

// Unmap then close the file
void writeFile_close()

#ifdef __MY_EXTERN_C__
}
#endif

#endif
