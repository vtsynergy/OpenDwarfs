/*
 ** CRC Kernel code
 **
 ** This code computes a 32-bit ethernet CRC using the "Slice-by-8" Algorithm published by Intel
 */

/* [FIX-OCL-03] was "../combinational-logic/crc/inc/eth_crc32_lut.h", which only resolved inside the old autotools build tree;
 * the header is now staged next to the kernel by the Makefile. */
#include "eth_crc32_lut.h"

__kernel void crc32_slice8(	__global const uint* restrict data, 
		uint length_bytes, 
		const uint length_ints,
		__global uint* restrict res,
		const uint num_items) /* [FIX-OCL-34] number of valid work-items (pages); global size is padded to a multiple of DWARF_LS */
{
	__private uint crc;
	__private uchar* currentChar;
	__private uint one,two;
	__private size_t i,j,gid;

	crc = 0xFFFFFFFF;
	gid = get_global_id(0);
	if (gid >= num_items) return; /* [FIX-OCL-34] guard the padded work-items */
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
		/* [FIX-OCL-03] was (unsigned char*): converts a private pointer to the generic address space, rejected by
		 * NVIDIA's OpenCL compiler ("changes address space of pointer") */
		currentChar = (__private uchar*) &one;
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
