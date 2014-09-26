#ifndef __PGI_SHIM_H
#define __PGI_SHIM_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  PGI_DEV_CPU=0,
  PGI_DEV_CUDA,
  PGI_DEV_OPENCL,
  PGI_DEV_NOT_HOST
} pgi_device_t;

uint32_t pgi_get_level();
uint32_t pgi_get_thread_num();
uint32_t pgi_get_num_threads();
uint32_t pgi_get_num_available_threads();

uint32_t pgi_get_device_num(pgi_device_t type);
uint32_t pgi_get_num_devices(pgi_device_t type);

void pgi_set_device(uint32_t id, pgi_device_t type);
void pgi_set_device_type(pgi_device_t type);
void pgi_set_device_by_id(uint32_t gpu_id_, pgi_device_t type);

void * pgi_get_opencl_queue(long async);
void * pgi_get_current_opencl_context();
void * pgi_get_current_opencl_device();

void pgi_sync(pgi_device_t type);

uint64_t pgi_get_memory();
uint64_t pgi_get_free_memory();

void pgi_run_parallel(int num_threads, void(*)(int, void*), void*);
void pgi_barrier();
void pgi_run_acc_region_to_init();

void* pgi_malloc(uint64_t size);
void pgi_free(void * buffer);
void pgi_memcpy_to(void *devptr, void* hostptr, uint64_t bytes);
void pgi_memcpy_from(void *hostptr, void* devptr, uint64_t bytes);

#ifdef __cplusplus
}
#endif

#endif /* __PGI_SHIM_H */

