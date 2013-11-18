/*
 ** CRC Kernel code
 **
 ** This code computes a 32-bit ethernet CRC using the "Slice-by-8" Algorithm published by Intel
 */

#include "../combinational-logic/crc/inc/eth_crc32_lut.h"

__kernel void crc32_slice8(	__global const uint* restrict data, 
		uint length_bytes, 
		const uint length_ints,
		__global uint* restrict res)
{
	__private uint crc;
	__private uchar* currentChar;
	__private uint one,two;
	__private size_t i,j,gid;

	crc = 0xFFFFFFFF;
	gid = get_global_id(0);
	i = gid * length_ints;

	while (length_bytes >= 8) // process eight bytes at once
	{
		one = data[i++] ^ crc;
		two = data[i++];
		crc = crc32Lookup[7][ one      & 0xFF] ^
			crc32Lookup[6][(one>> 8) & 0xFF] ^
			crc32Lookup[5][(one>>16) & 0xFF] ^
			crc32Lookup[4][ one>>24        ] ^
			crc32Lookup[3][ two      & 0xFF] ^
			crc32Lookup[2][(two>> 8) & 0xFF] ^
			crc32Lookup[1][(two>>16) & 0xFF] ^
			crc32Lookup[0][ two>>24        ];
		length_bytes -= 8;
	}

	while(length_bytes) // remaining 1 to 7 bytes
	{
		one = data[i++];
		currentChar = (unsigned char*) &one;
		j=0;
		while (length_bytes && j < 4) 
		{
			length_bytes = length_bytes - 1;
			crc = (crc >> 8) ^ crc32Lookup[0][(crc & 0xFF) ^ currentChar[j]];
			j = j + 1;
		}
	}

	res[gid] = ~crc;
}
