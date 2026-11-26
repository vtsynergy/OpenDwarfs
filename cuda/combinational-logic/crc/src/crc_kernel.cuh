/* CUDA translation of src/crc_kernel.cl ("Slice-by-8"): one page per thread, DWARF_LS threads/block, guard gid < num_items.
 * The lookup table of ../inc/eth_crc32_lut.h is included a second time as a __constant__ device array. */
#define OPENCL
#define __constant __constant__
#define crc32Lookup crc32Lookup_dev
#include "../inc/eth_crc32_lut.h"
#undef crc32Lookup
#undef __constant
#undef OPENCL

__global__ void crc32_slice8(	const unsigned int* __restrict__ data,
		unsigned int length_bytes,
		const unsigned int length_ints,
		unsigned int* __restrict__ res,
		const unsigned int num_items)
{
	unsigned int crc;
	unsigned char* currentChar;
	unsigned int one,two;
	size_t i,j,gid;

	crc = 0xFFFFFFFF;
	gid = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
	if (gid >= num_items) return;
	i = gid * length_ints;

	while (length_bytes >= 8) // process eight bytes at once
	{
		one = data[i++] ^ crc;
		two = data[i++];
		crc = crc32Lookup_dev[7][ one      & 0xFF] ^
			crc32Lookup_dev[6][(one>> 8) & 0xFF] ^
			crc32Lookup_dev[5][(one>>16) & 0xFF] ^
			crc32Lookup_dev[4][ one>>24        ] ^
			crc32Lookup_dev[3][ two      & 0xFF] ^
			crc32Lookup_dev[2][(two>> 8) & 0xFF] ^
			crc32Lookup_dev[1][(two>>16) & 0xFF] ^
			crc32Lookup_dev[0][ two>>24        ];
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
			crc = (crc >> 8) ^ crc32Lookup_dev[0][(crc & 0xFF) ^ currentChar[j]];
			j = j + 1;
		}
	}

	res[gid] = ~crc;
}
