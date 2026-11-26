/* CUDA-tree copy of the allocation helpers of opencl/include/common_args.c (non-FPGA paths). */
#include "common_args.h"

void check(int b, const char* msg)
{
	if(!b) { fprintf(stderr,"error: %s\n\n",msg); exit(-1); }
}
#define DWARF_NEW_ARRAY(name, type) \
void* name(const size_t N, const char* error_msg) { void* ptr = malloc(N * sizeof(type)); check(ptr != NULL, error_msg); return ptr; }
DWARF_NEW_ARRAY(char_new_array, char)
DWARF_NEW_ARRAY(int_new_array, int)
DWARF_NEW_ARRAY(long_new_array, long)
DWARF_NEW_ARRAY(float_new_array, float)
void* float_array_realloc(void* ptr, const size_t N, const char* error_msg)
{
	ptr = realloc(ptr, N * sizeof(float));
	check(ptr != NULL, error_msg);
	return ptr;
}
