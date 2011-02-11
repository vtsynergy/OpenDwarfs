#include "blast.h"

// Reads a vbyte into a 64-bit integer
unsigned char* vbyte_get64vbyte(unsigned char* cmem, uint8* n)
{
	int count = 0;

	*(n) = 0;
	while (((unsigned char)(*(cmem)) >= 128)) {
		*(n) |= (uint8)(*(cmem) & 0x7f) << count;
		count += 7;
		(cmem)++;
	}
	*(n) |= (uint8)(*(cmem)) << count;
	(cmem)++;

	return cmem;
}
