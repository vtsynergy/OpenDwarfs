/* vim: set sts=4 sw=4 expandtab:*/
/*
 * =====================================================================================
 *
 *       Filename:  splitter.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/31/11 15:24:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tom Scogland (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef __CTSAR_H
#define __CTSAR_H

#include <sched.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>
#include <errno.h>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif //_GNU_SOURCE
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <sched.h>
#include <utmpx.h>

#include "ctsar_dbg.h"

#ifdef __cplusplus
class ctsar;
#else
typedef void ctsar;
#endif
typedef void * ctsar_mem_handle;

typedef enum ctsar_dev_type{//bitfield components
    CTSAR_DEV_CPU=1,
    CTSAR_DEV_GPU=2
}ctsar_dev_type;


typedef enum ctsar_mem{
    CTSAR_MEM_NONE=0,
    CTSAR_MEM_PERSIST=1,
    CTSAR_MEM_PARTIAL=1<<1,//allow partial data movement based on iterator value
    CTSAR_MEM_INPUT=1<<2,  //copy in to devices
    CTSAR_MEM_OUTPUT=1<<3, //copy out from devices
    CTSAR_MEM_COLUMN=1<<4, //allow partial movement by columns                      
//pending release, these are just stubs for now
    CTSAR_MEM_PACKED=1<<5,  //attempt to pack data into smallest representation possible, rather than allocating complete region
    CTSAR_MEM_RESHAPE=1<<6  //attempt to reshape 2d segments into a better alignment
}ctsar_mem;



typedef struct ctsar_region{
    void * acc_data;
    void * host_data;

    size_t item_size;
    size_t rows;
    size_t cols;
    size_t device_pitch;
    size_t pad_x;
    size_t pad_y;
    int present;
    ctsar_mem flags;
}ctsar_region;

typedef void(* reduc_function)(void*,void*,void*);

typedef struct ctsar_reduction{
    void * reduc_var;//for reduction handling
    reduc_function reduc;
    size_t item_size;
    void * identity;
    void * host_data;
}ctsar_reduction;

typedef struct ctsar_device{
    ctsar_dev_type type;//type of device, for ratio calculations

    int start;//beginning of iteration space to execute in the next pass
    int end;//end of iteration space to execute in the next pass

    struct timeval time_start;
    struct timeval time_end;

    int active;//whether this device is to be used

    double pass_time;
    double time;
    double avg_tpi;//for dynamic chunk
    double tpi;

    int requested_iterations;//for dynamic chunk
    int iterations;
    int pass_iterations;
    int total_inner;//for dynamic chunk
    int inner;

    int group;

    size_t num_regions;
    ctsar_region regions[500];
    size_t free,total;
}ctsar_device;

typedef struct ctsar_dev_group{
    int members[50];//not safe, but expedient
    int size;
    double avg_tpi;//stores maximum tpi now, better scheduling results
    double part;
    int iterations;
    int fall_count;
    ctsar_dev_type type;
}ctsar_dev_group;

typedef enum ctsar_type{
    /* ctsar_NONE_CPU=0,
    ctsar_NONE_GPU=1, */
    CTSAR_NONE=0,
    CTSAR_STATIC=2,
    CTSAR_DYNAMIC=3,
    CTSAR_DEEP=4,
    CTSAR_START=5,
    CTSAR_CHUNK=6,
    CTSAR_CHUNK_STATIC=7,
    CTSAR_CHUNK_DYNAMIC=8,
    CTSAR_HYBRID_CHUNK=9//,
    /* ctsar_CHUNK_NOBAR=9,
    ctsar_CHUNK_STATIC_NOBAR=10,
    ctsar_CHUNK_DYNAMIC_NOBAR=11 */
}ctsar_type;


#ifdef __cplusplus
extern "C" {
#endif
void ctsar_next(ctsar* c, size_t size);
/* splitter * ctsar_next_chunk(int size, ctsar_dev_type dev_type); */
void ctsar_pre_init();
void ctsar_init( ctsar **s,
        int size,//total number of iterations to be split
        ctsar_type in_type,
        ctsar_dev_type devices,//bitfield, or the ENUM values to use >1
        double * in_ratio,//0-1 % to run on CPU pointer to allow for null
        int * in_depth,//also for null, optional arguments
        int * chunk_size//also for null, optional arguments
        );

void ctsar_end(ctsar*s);
void ctsar_start(ctsar*s);

int ctsar_get_type(ctsar*s);

int ctsar_loop(ctsar*s);
int ctsar_nested();

void assign_iterations(ctsar * s, int num_devs, int * devs, int * iters);

void mycg(int gts, int gte, double * ca, int * crs, int * cci, double * data, double * p);
void ctsar_retarget_mem(ctsar * c, void * oldp, void * newp);


void  * ctsar_reg_reduc(ctsar * c, void * ptr, void * identity, size_t item_size, reduc_function reduc); 
void  * ctsar_reg_mem(ctsar * c, void * ptr, size_t item_size, size_t count, 
                      ctsar_mem flags);
void *  ctsar_reg_mem_2d(ctsar * c, void * ptr, size_t item_size, size_t rows, size_t cols, 
                      ctsar_mem flags, size_t pad_x, size_t pad_y, size_t *pitch_output);
void    ctsar_unreg_mem(ctsar * c, void * ptr);
void    ctsar_swap_mem(ctsar * c, void * h1, void * h2);

//helpers
void **ctsar_alloc_2d(size_t dts, size_t sa, size_t sb);
void ***ctsar_alloc_3d(size_t dts, size_t sa, size_t sb, size_t sc);
void * ctsar_calloc_test(size_t size, size_t count);
void cMemGetInfo(size_t *free, size_t * total);

#define CSTART(X,Y) ctsar_get_start(X,Y)
#define CEND(X,Y) ctsar_get_end(X,Y)

uint64_t ctsar_get_start(ctsar *c, int tid);
uint64_t ctsar_get_end(ctsar *c, int tid);


#ifdef __cplusplus
}
#endif

#endif // __CTSAR_H
