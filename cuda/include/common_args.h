/* CUDA-tree replacement for opencl/include/common_args.h: only the host helpers the shared plain-C sources use (no OpenCL). */
#ifndef __COMMON_ARGS_H__
#define __COMMON_ARGS_H__
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#define MINIMUM(i,j) ((i)<(j) ? (i) : (j))
#ifdef __cplusplus
extern "C" {
#endif
extern void check(int b, const char* msg);
extern void* char_new_array(const size_t N, const char* error_msg);
extern void* int_new_array(const size_t N, const char* error_msg);
extern void* long_new_array(const size_t N, const char* error_msg);
extern void* float_new_array(const size_t N, const char* error_msg);
extern void* float_array_realloc(void* ptr, const size_t N, const char* error_msg);
#ifdef __cplusplus
}
#endif
#endif
