#include "common_args.h"
void check(int b, const char* msg) { if (!b) { fprintf(stderr, "error: %s\n\n", msg); exit(-1); } }
static void* checked(void* p, const char* msg) { check(p != NULL, msg); return p; }
void* char_new_array(const size_t N, const char* m)  { return checked(malloc(N * sizeof(char)), m); }
void* int_new_array(const size_t N, const char* m)   { return checked(malloc(N * sizeof(int)), m); }
void* long_new_array(const size_t N, const char* m)  { return checked(malloc(N * sizeof(long)), m); }
void* float_new_array(const size_t N, const char* m) { return checked(malloc(N * sizeof(float)), m); }
void* float_array_realloc(void* p, const size_t N, const char* m) { return checked(realloc(p, N * sizeof(float)), m); }
