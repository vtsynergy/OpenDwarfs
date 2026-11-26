/* Host-only subset of opencl/include/common_args.h, used by the host C sources copied from the OpenCL tree. */
#ifndef __COMMON_ARGS_H__
#define __COMMON_ARGS_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#ifdef __cplusplus
extern "C" {
#endif
#define MINIMUM(i,j) ((i)<(j) ? (i) : (j))
void check(int b, const char* msg);
void* char_new_array(const size_t N, const char* error_msg);
void* int_new_array(const size_t N, const char* error_msg);
void* long_new_array(const size_t N, const char* error_msg);
void* float_new_array(const size_t N, const char* error_msg);
void* float_array_realloc(void* ptr, const size_t N, const char* error_msg);
#ifdef __cplusplus
}
#endif
#endif
