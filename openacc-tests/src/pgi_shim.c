#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <accel.h>
#include <sched.h>
#include <errno.h>
#include <unistd.h>

#include "pgi_shim.h"

// forward decs for PGI opencl accessors
void* acc_get_opencl_queue(long async);
void* acc_get_current_opencl_context();
void* acc_get_current_opencl_device();

uint32_t pgi_get_level() { return omp_get_level(); }

uint32_t pgi_get_thread_num() { return omp_get_thread_num(); }

uint32_t pgi_get_num_threads() { return omp_get_num_threads(); }

uint32_t pgi_get_num_available_threads() { return omp_get_max_threads(); }

void* pgi_get_opencl_queue(long async) { return acc_get_opencl_queue(async); }
void* pgi_get_current_opencl_context() {
  return acc_get_current_opencl_context();
}
void* pgi_get_current_opencl_device() {
  return acc_get_current_opencl_device();
}

static acc_device_t translate_device_type(pgi_device_t type) {
  switch (type) {
    case PGI_DEV_CPU:
      return acc_device_host;
    case PGI_DEV_CUDA:
      return acc_device_nvidia;
    case PGI_DEV_OPENCL:
      return acc_device_radeon;
    case PGI_DEV_NOT_HOST:
      return acc_device_not_host;
    default:
      fprintf(stderr, "unknown device type (%d) passed to PGI shim, dying horribly", type);
      abort();
  }
}

static const char* translate_device_type_to_string(pgi_device_t type) {
  switch (type) {
    case PGI_DEV_CPU:
      return "CPU";
    case PGI_DEV_CUDA:
      return "CUDA";
    case PGI_DEV_OPENCL:
      return "radeon";
    case PGI_DEV_NOT_HOST:
      return "not host";
    default:
      return "unknown device";
  }
}

uint32_t pgi_get_device_num(pgi_device_t type) {
  return acc_get_device_num(translate_device_type(type));
}

uint32_t pgi_get_num_devices(pgi_device_t type) {
  return acc_get_num_devices(translate_device_type(type));
}

void pgi_set_device_type(pgi_device_t type) {
  acc_set_device_type(translate_device_type(type));
}

void pgi_set_device(uint32_t id, pgi_device_t type) {
  pgi_set_device_type(type);
  fprintf(stderr, "device type set to %d, num %u\n", acc_get_device_type(), id);
  int devid = acc_get_deviceid(id, translate_device_type(type));
  acc_set_deviceid(devid);
  pgi_run_acc_region_to_init();
}

void pgi_set_device_by_id(uint32_t gpu_id_, pgi_device_t type) {
  int deviceid = acc_get_deviceid(gpu_id_, translate_device_type(type));
  acc_set_deviceid(deviceid);
}

void pgi_sync(pgi_device_t type) { acc_sync(translate_device_type(type)); }

uint64_t pgi_get_memory() { return acc_get_memory(); }

uint64_t pgi_get_free_memory() { return acc_get_free_memory(); }

static int atomic_counter = 0;
void pgi_run_parallel(int num_threads, void (*fun)(int, void*), void* arg) {
  int siblings = pgi_get_num_threads();
// TODO: fix this for deeper nesting levels, potential race condition here
#pragma omp master
  { atomic_counter = 0; }
  int my_ticket = pgi_get_thread_num();

#pragma omp barrier
  while (my_ticket != atomic_counter) {
    sched_yield();
  }
#pragma omp parallel num_threads(num_threads)
  {
#pragma omp master
    {
#pragma omp atomic
      atomic_counter++;
      while (atomic_counter < siblings) {
        sched_yield();
      }
    }
    fun(pgi_get_thread_num(), arg);
  }
}

void pgi_barrier() {
#pragma omp barrier
}

#define MOCK_LOOP_SIZE 100
void pgi_run_acc_region_to_init() {
  int stuff[MOCK_LOOP_SIZE] = {0};
  int * stuff_d = (int*)pgi_malloc(sizeof(stuff[0])*MOCK_LOOP_SIZE);
  pgi_memcpy_to(stuff_d, stuff, sizeof(stuff[0])*MOCK_LOOP_SIZE);
#pragma acc kernels loop independent deviceptr(stuff_d)
  for (int i = 0; i < MOCK_LOOP_SIZE; ++i) {
    stuff_d[i] = i;
  }
  pgi_memcpy_from(stuff, stuff_d, sizeof(stuff[0])*MOCK_LOOP_SIZE);
  pgi_free(stuff_d);
}

void* pgi_malloc(uint64_t size) {
  fprintf(stderr, "PGI device type %d id %d allocating size %zu, %lu from %zu free of %zu total\n", acc_get_device_type(), acc_get_device_num(acc_get_device_type()), size, (__PGI_ULLONG)size, pgi_get_free_memory(), pgi_get_memory());
  void* ret = acc_malloc(size);
  if (ret == NULL) {
    fprintf(stderr, "Device allocation failed, requested size %zu\n", size);
    exit(ENOMEM);
  } else {
    if (ret == (void*)0x80) {
      fprintf(stderr, "PGI allocation invalid, retrying\n");
      pgi_free(ret);
      ret = acc_malloc((__PGI_ULLONG)size);
      if (ret == NULL) {
        fprintf(stderr, "Device allocation failed, requested size %zu\n", size);
        exit(ENOMEM);
      } else if (ret == (void*)0x80) {
        fprintf(stderr, "Device allocation failed, bad address\n");
        exit(ENOMEM);
      }
    }
    /* fprintf(stderr, */
    /*         "PGI allocation successful, size %zu at location %p\n", */
    /*         size, */
    /*         ret); */
  }
  return ret;
}
void pgi_memcpy_to(void* devptr, void* hostptr, uint64_t bytes) {
  /* fprintf(stderr, */
  /*         "PGI memcpy to: dev dst:%p host src:%p bytes:%zu\n", */
  /*         devptr, */
  /*         hostptr, */
  /*         bytes); */
  acc_memcpy_to_device(devptr, hostptr, bytes);
}
void pgi_memcpy_from(void* hostptr, void* devptr, uint64_t bytes) {
  /* fprintf(stderr, */
  /*         "PGI memcpy from: host dst:%p dev src:%p bytes:%zu\n", */
  /*         hostptr, */
  /*         devptr, */
  /*         bytes); */
  acc_memcpy_from_device(hostptr, devptr, bytes);
}

void pgi_free(void* buffer) { acc_free(buffer); }
