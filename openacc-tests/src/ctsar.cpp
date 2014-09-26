/* vim: set sts=4 sw=4 expandtab:*/
/* 
 *  Copyright 2010 by Virginia Polytechnic Institute and State
 *  University. All rights reserved. Virginia Polytechnic Institute and
 *  State University (Virginia Tech) owns the software and its
 *  associated documentation.
 *
 *  Contact: Tom Scogland <tom.scogland@gmail.com>
 */

#include <stdio.h>
#include "ctsar_impl.hpp"
#include "pgi_shim.h"
#include <sched.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))
#define max(x1,x2) ((x1) > (x2) ? (x1):(x2))


#define DEFAULT_CHUNK 1000
#define DEFAULT_DEPTH 10


using namespace std;

//global flags, all apply across regions
int sdebug=-1;
int allowed_devs;
int stop_slow=1;
int nested=0;//whether to group all devices of a type and treat them as one, how the original version worked

// once_flag topology_once;
// hwloc_topology_t topology;


// factoring out boost tss
inline unsigned int get_tid(){
    // return *my_tid;
    return pgi_get_thread_num();
}
inline void set_tid(unsigned int new_tid){
    // *my_tid = new_tid;
}
inline void init_tid(){
    // if(!my_tid.get()){
    //     my_tid.reset(new unsigned int(pgi_get_thread_num()));
    // }
}

static struct timeval start_time;

void cleanup() {
    fflush(stdout);
    fflush(stderr);
    sleep(1);
}

static struct timeval ctsar_time_start;
static struct timeval ctsar_time_end;

static void init_time() {
    gettimeofday(&start_time, NULL);
}

static uint64_t get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (uint64_t) (t.tv_sec - start_time.tv_sec) * 1000000
        + (t.tv_usec - start_time.tv_usec);
}

uint64_t calc_time(struct timeval start, struct timeval t) {
    return (uint64_t) (t.tv_sec - start.tv_sec) * 1000000
        + (t.tv_usec - start.tv_usec);
}


ctsar::~ctsar(){
}
ctsar::ctsar(){
}


extern "C" {

void ctsar_print_messages(ctsar *s){
    // s->all_barr->wait();
    pthread_barrier_wait(&s->all_barr);
    if((get_tid()) == 0){
        s->outer_iteration++;
        if(sdebug > DBG_TEST)
            fprintf(stderr, "ctsar:\tpass:%zd\tsched_time:%ld\n", s->outer_iteration, s->ctsar_next_time);
        s->fully_assigned=0;
        s->assigned = 0;
    }
    // s->all_barr->wait();
        pthread_barrier_wait(&s->all_barr);
    ctsar_device *d = &s->old_devices[get_tid()];
    if(sdebug >= 0){
        fprintf(stderr, "ctsar:\tpass:%lu\tinner:%d\tdevice:%d\ttype:%s\ti:%d\tt:%e\ttpi:%e\n", s->outer_iteration, d->inner, get_tid(), d->type == CTSAR_DEV_CPU ? "CPU" : "GPU", d->pass_iterations, d->time, d->time/d->pass_iterations);
        // s->all_barr->wait();
        pthread_barrier_wait(&s->all_barr);
        fprintf(stderr, "OMP_CTSAR_RATIO_%d=%lf\n", get_tid(), d->time/d->pass_iterations);
    }
    s->old_devices[get_tid()].time = 0;
    s->old_devices[get_tid()].inner = 0;
        /* for(i=0; i<s->num_devices; ++i){ */
            s->old_devices[get_tid()].pass_iterations = 0;
        /* } */
    // s->all_barr->wait();
        pthread_barrier_wait(&s->all_barr);
}

void set_gpu(int itid){
    if(pgi_get_num_devices(PGI_DEV_NOT_HOST) > 0){
        uint32_t request = itid % pgi_get_num_devices(PGI_DEV_NOT_HOST);
        uint32_t num_cuda_devices = pgi_get_num_devices(PGI_DEV_CUDA);
        uint32_t num_radeon_devices = pgi_get_num_devices(PGI_DEV_OPENCL);
        if(request < num_cuda_devices) {
            pgi_set_device(request, PGI_DEV_CUDA);
        }else if(request < num_cuda_devices+num_radeon_devices){
            pgi_set_device(request - num_cuda_devices, PGI_DEV_OPENCL);
        }else{
            fprintf(stderr, "unknown non-host device type found, aborting\n");
            abort();
        }
        /*acc_get_current_cuda_context();*/
        int dev = pgi_get_device_num(PGI_DEV_NOT_HOST);
        DBG(DBG_DEBUG,"[%d] selected device number %d of %d, got %d/%d\n", itid, request, pgi_get_num_devices(PGI_DEV_NOT_HOST),pgi_get_device_num(PGI_DEV_NOT_HOST), dev);
    }
}

void reg_gpu(ctsar * s){
    static int initialized[50] = {0};//TODO: not safe for systems with many devices
    if(!initialized[get_tid()]){
        if(s->old_devices[get_tid()].type == CTSAR_DEV_GPU){
            set_gpu(get_tid());
            ctsar_device *dev = &s->old_devices[get_tid()];
            dev->total = pgi_get_memory();
            dev->free = pgi_get_free_memory();
        }else{
            pgi_set_device_type(PGI_DEV_CPU);
        }
        initialized[get_tid()]++;
    }
}

void ctsar_pre_init(){
    // fprintf(stderr, "There are %d core objects in the topology\n",hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE));

    // this is here to ensure that a single thread checks this first, having
    // all threads do it can cause a deadlock in PGI's init routines
    int unnecessary_nd = pgi_get_num_devices(PGI_DEV_NOT_HOST);
    pgi_run_acc_region_to_init();

    pgi_run_parallel(pgi_get_num_available_threads(), [](int omp_tid, void * unnecessary){
        int nd = pgi_get_num_devices(PGI_DEV_NOT_HOST);
        init_tid();
        // fprintf(stderr,"My get_tid() is %d\n", get_tid());
        char * env = NULL;
        if(!(((env = getenv("CTSAR_NOBIND")) != NULL) && atoi(env))){
        //hwloc and libnuma from PGI's openacc do not play nice
            // pid_t ptid = syscall(SYS_gettid);
            // {
            //     hwloc_obj_t core = hwloc_get_obj_by_type(topology,HWLOC_OBJ_CORE,get_tid());
            //     if(core == NULL){
            //         fprintf(stderr,"Not good... %d\n", get_tid());
            //     }else{
            //         fprintf(stderr,"Success..? %d\n", get_tid());
            //     hwloc_cpuset_t set = hwloc_bitmap_dup(core->cpuset);
            //     fprintf(stderr,"binding thread: %d\n", get_tid());
            //     hwloc_set_cpubind(topology, set, HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT);
            //     hwloc_bitmap_free(set);
            //     }
            // }
            pid_t ptid = 0;
            cpu_set_t cpus;

            CPU_ZERO(&cpus);
            CPU_SET(get_tid(),&cpus);
            // fprintf(stderr,"omp tid:%d, num threads\n",get_tid(), pgi_get_num_threads());
            sched_setaffinity(0, sizeof(cpu_set_t), &cpus);
            sched_getaffinity(0, sizeof(cpu_set_t), &cpus);
            for(int i=0; i < pgi_get_num_threads(); i++){
                if(CPU_ISSET(i, &cpus))
                    fprintf(stderr, "tid:%d, has cpu:%d set\n", get_tid(), i);
            }
        }
        if(get_tid() < nd){
            set_gpu(get_tid());
            //acc_init(acc_device_nvidia);
        }else{
            //acc_init(acc_device_host);
        }
    }, nullptr);
}

void ctsar_cpu_barrier(ctsar *s){
    if(ctsar_get_type(s)==CTSAR_DEV_CPU)
    {
        // s->b_barr->wait();
        pthread_barrier_wait(&s->b_barr);
        //pthread_barrier_wait(&s->barr);
    }
}

int ctsar_get_type(ctsar*s)
{
    return s->old_devices[get_tid()].type;
}

int update_device(ctsar*s, ctsar_mem flags){
    int sync=0;
    int i;
    ctsar_device * dev = &s->old_devices[get_tid()];
    reg_gpu(s);
    if(dev->type != CTSAR_DEV_CPU){
        for(i=0; i < dev->num_regions; i++){
            ctsar_region * reg = &dev->regions[i];
            if(  (reg->acc_data == NULL) //blank region
                    ||( ! (reg->flags & CTSAR_MEM_INPUT))//not input
                    ||((reg->flags & flags) != flags)//does not match flags
                    ||( reg->present && (!( reg->flags & CTSAR_MEM_PARTIAL))))//already copied
                continue;


            size_t size = reg->item_size * reg->rows * reg->cols;
            size_t host_offset = 0;
            size_t device_offset = 0;
            size_t iterations = (dev->end - dev->start);
            size_t row_range = reg->rows;

            size_t row_width = reg->item_size * reg->cols;
            if(reg->flags & CTSAR_MEM_PARTIAL){
                size_t py_s = (((signed)dev->start - (signed)reg->pad_y) < 0) ? 0 : reg->pad_y;
                size_t py_e = ((dev->end + reg->pad_y) > reg->rows) ? 0 : reg->pad_y;
                host_offset = (dev->start-py_s) * reg->item_size * reg->cols;
                device_offset = host_offset;
                /* fprintf(stderr,"%d: py_s:%d py_e:%d size:%lu copy size:%lu copy offset:%lu start:%d end:%d\n", get_tid(), py_s, py_e, size, new_size, offset, dev->start, dev->end); */
                size = (iterations + (py_s+py_e)) * row_width;
                row_range = iterations;
            }

            if(reg->flags & CTSAR_MEM_COLUMN){//copy data to GPU partial and column-wise
                size_t csize   = reg->item_size * iterations;
                size_t coffset =  reg->item_size * dev->start;
                DBG(DBG_ALWAYS, "copying %lu bytes to/from GPU %d as a column\n", reg->rows * csize, get_tid());
                for(int j=0; j<reg->rows; ++j){
                    if(reg->flags & CTSAR_MEM_PARTIAL//skip column rows which are copied by partial
                            &&(coffset > device_offset && coffset+csize < (device_offset+size))) {
                        /* DBG(DBG_ALWAYS, "skipping row %d of %d to/from GPU %d as part of a column\n", j, reg->rows, get_tid()); */
                    }else{
                        /* DBG(DBG_ALWAYS, "copying row %d of %d to/from GPU %d as part of a column\n", j, reg->rows, get_tid()); */
                        pgi_memcpy_to(((char*)reg->acc_data) + coffset, ((char*)reg->host_data) + coffset, csize);
                    }
                    coffset += (reg->cols * reg->item_size);
                }
                reg->present=1;
                if(! (reg->flags & CTSAR_MEM_PARTIAL) ){//do not do copy below undless the partial flag is set
                    continue;
                }
            }

            int dev;
            dev = pgi_get_device_num(PGI_DEV_NOT_HOST);
            DBG(DBG_DEBUG, "copying %lu bytes to GPU %d, addr=%p, host offset=%lu, device offset=%lu\n", size, get_tid(), reg->acc_data, host_offset, device_offset);
            pgi_memcpy_to(((char*)reg->acc_data) + device_offset, ((char*)reg->host_data) + host_offset, size);

            sync=1;
            reg->present=1;
        }
    }
    return sync;
}

int update_host(ctsar*s, ctsar_mem flags){
    int i;
    ctsar_device * dev = s->old_devices+get_tid();
    int sync=0;

    if(dev->type != CTSAR_DEV_CPU){
        for(i=0; i < dev->num_regions; i++){
            ctsar_region * reg = &dev->regions[i];
            if(reg->acc_data == NULL)
                continue;
            if( ! (reg->flags & CTSAR_MEM_OUTPUT))//not output
                continue;
            if((reg->flags & flags) != flags)//does not match flags
                continue;

            size_t size = reg->item_size * reg->rows * reg->cols;
            size_t row_range = reg->rows;
            size_t offset = 0;
            size_t device_offset = 0;
            if(reg->flags & CTSAR_MEM_PARTIAL){
                row_range = dev->end - dev->start;
                offset = dev->start * reg->item_size * reg->cols;
                device_offset = offset;
                size = row_range * reg->item_size * reg->cols;
            }

            size_t row_width = reg->item_size * reg->cols;

            /* DBG(-1, "in end, found region to test, flags=%d\n",reg->flags); */
            DBG(DBG_DEBUG, "[%u] copying from gpu, offset=%lu, size=%lu\n",get_tid(),offset,size);
            pgi_memcpy_from(((char*)reg->host_data)+offset, ((char*)reg->acc_data) + device_offset, size);
            sync++;
        }
    }
    return sync;
}

void ctsar_start(ctsar*s){
    ctsar_device * dev = &s->old_devices[get_tid()];

    if(s->old_devices[get_tid()].type == CTSAR_DEV_GPU){
        update_device(s,CTSAR_MEM_PERSIST);
        pgi_sync(PGI_DEV_NOT_HOST);
    }

    gettimeofday(&(dev->time_start), NULL);

    if(s->old_devices[get_tid()].type == CTSAR_DEV_GPU){
        update_device(s,CTSAR_MEM_NONE);
        pgi_sync(PGI_DEV_NOT_HOST);
    }

    if(s->type >= CTSAR_CHUNK) {//time copies only for non-chunk schedulers
        /* pgi_sync(PGI_DEV_NOT_HOST); */
        gettimeofday(&(dev->time_start), NULL);
    }


    //printf("%d: starting loop it=%d\n", omp_get_thread_num(), s->iteration);
}


void ctsar_end(ctsar*s) {
    /* if(nested && s.old_devices[get_tid()].type == ctsar_DEV_CPU){
        ctsar_cpu_barrier();
    } */
    //printf("%d: ending loop it=%d\n", omp_get_thread_num(), s->iteration);
    uint64_t time;
    ctsar_device * dev = s->old_devices+get_tid();


    if(s->old_devices[get_tid()].type == CTSAR_DEV_GPU){
        pgi_sync(PGI_DEV_NOT_HOST);
        if(update_host(s, CTSAR_MEM_NONE) && s->type < 6)
            pgi_sync(PGI_DEV_NOT_HOST);
    }
    gettimeofday(&(dev->time_end), NULL);

    time = calc_time(dev->time_start,dev->time_end);

    dev->pass_time = (time/1000000.0);
    dev->time += (double)time;
    dev->pass_iterations += dev->iterations;
    /* if(dev->type == CTSAR_DEV_CPU)
        dev->pass_time *= 1.05; */

    if(dev->iterations > 0){
        /* dev->tpi = (((dev->pass_time / dev->iterations) + dev->tpi) / 2.0); */
        dev->tpi = (time/1000000.0) / dev->iterations;

        /* dev->avg_tpi = ((dev->avg_tpi * (dev->total_inner - 1)) + dev->tpi) / dev->total_inner; */
        if(dev->avg_tpi != 0){
            dev->avg_tpi = ((dev->avg_tpi * 5) + dev->tpi) / 6;
        }else{
            dev->avg_tpi = dev->tpi;
        }
    }
    /* if(get_tid()==0) */
    DBG(DBG_DEBUG,"el[%d]: it=%d, loop_time=%lu, tpi=%lf, tpie=%e, tpia=%lf\n", get_tid(), dev->iterations, time, dev->tpi, dev->tpi, dev->avg_tpi);
}		/* -----  end of function ctsar_cpu_end  ----- */

void ctsar_init(ctsar**c,
        int size,//total number of iterations to be split
        ctsar_type in_type,
        ctsar_dev_type devices,//bitfield, or the ENUM values to use >1
        double * in_ratio,//0-1 % to run on CPU pointer to allow for null
        int * in_depth,//also for null, optional arguments
        int * chunk_size//also for null, optional arguments
        )
{
    struct timeval init_start, init_end;

    gettimeofday(&init_start, NULL);
    char *env;
    int i;
    ctsar_dev_type dev = devices;
    int threads = pgi_get_num_available_threads();

    if(*c == NULL && pgi_get_thread_num() == 0){
    fprintf(stderr, "%p\n", *c);
        *c = new ctsar();
        ctsar *s = *c;
        memset(s, 0, sizeof(ctsar));
        if(s == NULL){
            fprintf(stderr, "could not allocate ctsar construct, failing\n");
            exit(1);
        }
        // Default values
        s->type = CTSAR_DYNAMIC;
        s->depth = 10;
        s->chunk = 100;

        {//arguments/environment
            {//set type
                if((env = getenv("OMP_CTSAR_TYPE")) != NULL){
                    s->type = (ctsar_type)atoi(env);
                }else{
                    s->type = CTSAR_DYNAMIC;
                }

            }
            {//set depth
                if((env = getenv("OMP_CTSAR_DEPTH")) != NULL){
                    s->depth = atoi(env);
                }else if(in_depth != NULL){
                    s->depth = *in_depth;
                }else{
                    s->depth = 10;
                }
            }
            {//set chunk size
                if((env = getenv("OMP_CTSAR_CHUNK")) != NULL){
                    s->chunk = (ctsar_type)atoi(env);
                }else if(chunk_size != NULL){
                    s->chunk = *chunk_size;
                }else{
                    s->chunk = 100; //generic default
                }

                s->chunk = size / s->chunk;
                if (s->chunk < 1)
                    s->chunk = 1;
            }
            {//set flags
                if((env = getenv("OMP_CTSAR_DBG")) != NULL)
                    sdebug = atoi(env);
                else
                    sdebug = 0;

                if((env = getenv("OMP_CTSAR_NEST")) != NULL)
                    nested = atoi(env);
                else
                    nested = 0;

                if((env = getenv("CTSAR_STOP_SLOW")) != NULL){
                    stop_slow = atoi(env);
                }else{
                    stop_slow = 1;
                }
            }
            {//set ratio
                if(in_ratio != NULL){
                    s->ratio = *in_ratio;
                }else{
                    if(pgi_get_num_devices(PGI_DEV_NOT_HOST) > 0){
                        // double gpucores = cgetcores(0);
                        double gpucores = 448; //FIXME: static number of cores in 2070 for now
                        /* fprintf(stderr, "gpucores=%lf\n", gpucores); */
                        /* exit(1); */
                        if(nested){
                            s->ratio = 1.0/(gpucores/ //cores in a C2050 GPU
                                    //cores left on CPU
                                    4.0//float width per core
                                    );//split problem 1 per sse lane to 1 per GPU core (448)
                        }else{
                            s->ratio = 1.0/(gpucores/ //cores in a C2050 GPU
                                    4.0//float width per core
                                    );//split problem 1 per sse lane to 1 per GPU core (448)
                        }
                        if(s->type >= CTSAR_CHUNK){
                            s->ratio *= 10; //reduce the initial difference between CPU and GPU
                        }
                    }else{
                        s->ratio = 1;
                    }
                }
                if((env = getenv("OMP_CTSAR_RATIO")) != NULL)
                    sscanf(env, "%lf", &s->ratio);
                if(s->ratio > 1 || s->ratio < 0){
                    fprintf(stderr,"ratio is outside the range of 0-1, failing\n");
                    cleanup();
                    exit(1);
                }
            }
            {//devices allowed
                if((env = getenv("OMP_CTSAR_DEVICES")) != NULL)
                    dev = (ctsar_dev_type)atoi(env);

                allowed_devs = dev;
                if( ! ((dev & CTSAR_DEV_CPU) || (dev & CTSAR_DEV_GPU)) ){
                    fprintf(stderr, "Some device must be enabled, exiting\n");
                    cleanup();
                    exit(1);
                }
            }
        }


        {//initialize ctsar struct
            int gpus = 0, cpus;
            int cur_group=0;
            cerr << "Initializing CTSAR" << endl;

            if(dev & CTSAR_DEV_GPU){
                int num_gpus = pgi_get_num_devices(PGI_DEV_NOT_HOST);
                if(num_gpus < 1){
                    fprintf(stderr, "no GPUs found, failing\n");
                    cleanup();
                    exit(1);
                }
                if(((env = getenv("OMP_CTSAR_MULTI_GPU")) != NULL)
                        || ((env = getenv("CTSAR_MULTI_GPU")) != NULL)){
                    if(atoi(env) > 0){
                        if(num_gpus < atoi(env)){
                            gpus = num_gpus;
                        }else{
                            gpus = atoi(env);
                        }
                    }
                }else{
                    gpus = pgi_get_num_devices(PGI_DEV_NOT_HOST);
                }
            }else{
                gpus = 0;
            }

            if(dev & CTSAR_DEV_CPU){
                cpus=threads - gpus;//we make an assumption that the number of threads equals the number of cores for now, also that we assign a core to do *nothing* but manage GPUs
                if(cpus == 0){
                    fprintf(stderr, "CPUs and GPUs selected, no cores left on the CPU, failing\n");
                    cleanup();
                    exit(1);
                }

            }else{
                cpus = 0;
            }


            s->num_devices = threads;
            s->num_dev_groups = gpus+cpus;
            {//validity checks, some types only allowed with certain configurations
                if(s->type == CTSAR_CHUNK_DYNAMIC && s->num_dev_groups < 2){
                    fprintf(stderr,"The chunk dynamic scheduler requires two or more devices, falling back on chunk static\n");
                    s->type = CTSAR_CHUNK_STATIC;
                }
            }

            //allocate cpu range arrays based on the max number of threads
            /* s->old_devices = (ctsar_device *)calloc(sizeof(ctsar_device),threads); */
            s->old_devices = (ctsar_device *)calloc(sizeof(ctsar_device),pgi_get_num_available_threads());
            s->dev_groups = (ctsar_dev_group *)calloc(sizeof(ctsar_dev_group),s->num_dev_groups);

            for(i=0; i<gpus; i++){
                s->old_devices[i].type = CTSAR_DEV_GPU;
                s->old_devices[i].avg_tpi = 0;
                s->old_devices[i].tpi = s->ratio;//derived from original equation in ratio setting above
                s->old_devices[i].active = 1;
                s->old_devices[i].group = cur_group;

                s->dev_groups[cur_group].iterations = 1;

                s->dev_groups[cur_group].members[0] = i;
                s->dev_groups[cur_group].size = 1;
                s->dev_groups[cur_group].type = CTSAR_DEV_GPU;
                cur_group++;
            }

            for(;i<threads;i++){
                s->old_devices[i].type = CTSAR_DEV_CPU;
                s->old_devices[i].avg_tpi = 0;
                s->old_devices[i].tpi = 1-s->ratio;
                s->old_devices[i].start = 0;
                s->old_devices[i].end = 0;

                if(i< cpus+gpus){
                    s->old_devices[i].active = 1;
                }else{
                    s->old_devices[i].active = 0;
                }
                /* printf("setting up dev %d, type %d\n", i, s->old_devices[i].type); */
            }
            char env_var[50] = {0};
            for(i=0; i<threads; i++){
                sprintf(env_var, "OMP_CTSAR_RATIO_%d", i);
                if((env = getenv(env_var)) != NULL){
                    s->old_devices[i].tpi = atof(env);
                }
            }

            if(dev & CTSAR_DEV_CPU){
                s->dev_groups[cur_group].type = CTSAR_DEV_CPU;
                if(nested){
                    s->dev_groups[cur_group].size = threads-gpus;
                    for(i=0; i<threads-gpus; i++){
                        s->dev_groups[cur_group].members[i] = gpus+i;
                        s->old_devices[gpus+i].group = cur_group;
                    }
                    cur_group++;
                }else{
                    for(i=0; i<cpus; i++){
                        s->dev_groups[cur_group].type = CTSAR_DEV_CPU;
                        s->dev_groups[cur_group].members[0] = gpus+i;
                        s->dev_groups[cur_group].size = 1;
                        s->old_devices[gpus+i].group = cur_group;
                        cur_group++;
                    }
                }
            }

            if(threads - gpus > 0){
                pthread_barrier_init(&s->b_barr, NULL, threads-gpus);
                // s->b_barr = make_shared<boost::barrier>(threads-gpus);
            }
            fprintf(stderr, "#####all_barr initialized to %d threads\n", threads);
            pthread_barrier_init(&s->all_barr, NULL, threads);
            // s->all_barr = make_shared<boost::barrier>(threads);
            //pthread_barrier_init(&s->barr, NULL, threads-gpus);
            DBG(0,"SPLIT: initial parameters, type=%d, ratio=%lf, depth=%d, nested=%d, threads=%d, gpus=%d\n", s->type, s->ratio, s->depth,nested,threads, gpus);
        }
        s->outer_iteration=0;
    }

    //ignoring repeat initialization
    (*c)->pass_iterations = size;
    (*c)->iteration=0;
    /* for(i=0; i<s->num_devices; ++i){
       s->old_devices[i].pass_iterations = 0;
       } */
    gettimeofday(&init_end, NULL);
    DBG(DBG_DEBUG, "ctsar_init took %lu\n", calc_time(init_start,init_end));
}

void calculate_multi_ratio(ctsar * s, int iterations, int assign){
    int i;
    if(s->num_dev_groups == 1){
        int dev = 0;
        assign_iterations(s,1, &dev, &iterations);
    }else{
        ctsar_dev_group * ref_group = NULL;
        ctsar_dev_group * group;
        double lowest_tpi = 5000000; 
        for(i=0; i < s->num_dev_groups; i++){
            group = &s->dev_groups[i];
            group->avg_tpi = 0.0;
            for(int j=0; j < group->size; j++){
                if(s->old_devices[group->members[j]].tpi > group->avg_tpi)
                    group->avg_tpi = s->old_devices[group->members[j]].tpi;//using max to lessen blocking time
            }
            if((ref_group == NULL || group->avg_tpi <= lowest_tpi) && group->type == CTSAR_DEV_CPU)
                ref_group = group;
            if(group->avg_tpi < lowest_tpi)
                lowest_tpi = group->avg_tpi;
        }

        if(allowed_devs & CTSAR_DEV_CPU &&  stop_slow){
            for(i=0; i < s->num_dev_groups; i++){
                group = &s->dev_groups[i];
                if(group->type != CTSAR_DEV_CPU){
                    if(group->avg_tpi > ref_group->avg_tpi*1.1){
                        DBG(DBG_DEBUG,"STOP_SLOW: active, checking, me=%e cpus=%e\n",group->avg_tpi, ref_group->avg_tpi);
                        group->fall_count++;

                        if(group->fall_count>1){
                            DBG(DBG_ALWAYS,"STOP_SLOW: converting group=%d to CPU threads\n",i);
                            ctsar_device * dev;
                            for(int j=0; j<group->size; j++){
                                dev = &s->old_devices[group->members[j]];
                                dev->tpi = ref_group->avg_tpi;
                                dev->type = CTSAR_DEV_CPU;
                            }
                            group->avg_tpi = ref_group->avg_tpi;
                            group->type = CTSAR_DEV_CPU;

                            int cnt=0;
                            for(int j=0; j<s->num_devices;j++){
                                if(s->old_devices[j].type == CTSAR_DEV_CPU){
                                    cnt++;
                                }
                            }
                            if(cnt > 0){
                                // s->b_barr = make_shared<boost::barrier>(cnt);
                                pthread_barrier_init(&s->b_barr, NULL, cnt);
                            }
                            //pthread_barrier_init(&s->barr, NULL, cnt);
                        }
                    }else{
                        group->fall_count=0;
                    }
                }
            }
        }

        //solve for new ratio
        DBG(DBG_DEBUG,"solving for new ratios\n");
        if(sdebug) init_time();

        double total_dist = 0;
        for(i=0; i<s->num_dev_groups; i++){
            group = &s->dev_groups[i];
            double inverse = 1 / (group->avg_tpi/(double)group->size);
            group->part = inverse * iterations;
            total_dist += group->part;
        }

        DBG(DBG_DEBUG,"ratio calc time:%lu\n", get_time());

        //pull out the results
        int devs[s->num_dev_groups];
        int its[s->num_dev_groups];

        int total=0;
        for(i=0; i<s->num_dev_groups; i++){
            group = &s->dev_groups[i];
            devs[i]=i;
            if((int)((group->part / total_dist)*iterations) > 0){
                its[i] = (group->part / total_dist)*iterations;
            }else{
                its[i] = 1;
            }

            // fprintf(stderr, "ESTIMATING[%d]: %lf %lf %lf %lf\n", i, vars[i], (group->part / total_dist), (group->part / total_dist) * iterations, vars[i] * iterations);

            DBG(DBG_DEBUG,"group[%d]: ratio=%lf, its=%d, pred. time=%lf\n", i, (group->part / total_dist), its[i], its[i] * s->old_devices[i].tpi);
            total+= its[i];
        }

        if(total < iterations){
            i=0;
            while(total < iterations){
                if(s->dev_groups[i%s->num_dev_groups].type == CTSAR_DEV_GPU || its[i%s->num_dev_groups] >= 1000 || pgi_get_num_devices(PGI_DEV_NOT_HOST) <= 0){
                    its[i%s->num_dev_groups ]++;
                    total++;
                }
                if(i > s->num_dev_groups){
                    its[i%s->num_dev_groups ]++;
                    total++;
                }
                i++;
            }
        } else if(total > iterations)
            for(i=0; i < total - iterations; i++)
                its[i%s->num_dev_groups]--;

        if (assign)
            assign_iterations(s,s->num_dev_groups, devs,its);
    }
}
int get_remaining(ctsar * s){
    return s->pass_iterations - s->assigned;
}

void assign_iterations(ctsar * s, int num_groups, int * groups, int * iters)
{
    int i;

    /* DBG(1,"ndev:%d\n", s->num_devices); */

    for(i=0; i<num_groups; i++){
        int j, ctr;
        ctsar_dev_group * group = &s->dev_groups[groups[i]];
        group->iterations = iters[i];
        for(j=0, ctr=0; j< group->size; j++, ctr++){
            ctsar_device * dev = &s->old_devices[group->members[j]];

            int to_add = iters[i]/group->size;
            /* printf("to_add: %d size:%d iters[i]: %d\n", to_add, s->dev_groups[i].size, iters[i]); */
            if(iters[i] % group->size > ctr)
                to_add++;

            /* if(ctr == s->dev_groups[groups[i]].size-1 && dev->end + to_add < ) */

            dev->start = s->assigned;
            dev->end = s->assigned + to_add;
            if(dev->end > s->pass_iterations){
                dev->end = s->pass_iterations;
            }
            /* if(sdebug > 2) fprintf(stderr,"assigning thread %d %d more iterations, END: = %d\n", j, dev->end - s->assigned, dev->end); */
            s->assigned = dev->end;
            dev->iterations = dev->end - dev->start;
            DBG(DBG_DEBUG,"ai[%d]: group=%d, group_iters=%d, assigned=%d, start=%d, end=%d, iters=%d\n",
                    group->members[j], i, iters[i], dev->end - dev->start,dev->start, dev->end,dev->iterations);
        }
    }
    /* if(sdebug > 1) fprintf(stderr,"END: low end %d, high start %d\n", s->assigned, high_start); */
    if(s->assigned == s->pass_iterations){
        DBG(DBG_DEBUG,"all iterations assigned\n");
        s->fully_assigned=1;
    }
}

void ctsar_next(ctsar * s, size_t size //valid only for non-deep types
        ){
    size_t it=size;
        // fprintf(stderr, "thread %d of %u entering ctsar_next\n", get_tid(), pgi_get_num_threads());


    if( get_tid() == 0){
        DBG(DBG_DEBUG,"calculating next pass, iteration=%zu\n",s->iteration);
    }

    if(s->old_devices[get_tid()].type == CTSAR_DEV_GPU){
        reg_gpu(s);
    }
    if(get_tid()==0)
    {
        gettimeofday(&ctsar_time_start, NULL);
    }
    if(s->type < CTSAR_CHUNK){//if using a barrier method, block on a barrier
        //ensures that all tpi and time values are up to date before solve
        pthread_barrier_wait(&s->all_barr);
        // s->all_barr->wait();
        // fprintf(stderr, "thread %d entering next barrier\n", get_tid());

        if(get_tid() == 0)
        {
            int i;
            it=size;

            switch(s->type){
                case CTSAR_STATIC:
                    // USE initial values, do not redefine
                    for(i=0; i<s->num_devices; i++){
                        char env_var[50], *env;
                        sprintf(env_var, "OMP_CTSAR_RATIO_%d", i);
                        if((env = getenv(env_var)) != NULL){
                            s->old_devices[i].tpi = atof(env);
                        }else if(s->old_devices[i].type == CTSAR_DEV_CPU){
                            s->old_devices[i].tpi = 1-s->ratio;
                        }else{
                            s->old_devices[i].tpi = s->ratio;
                        }
                    }
                    break;
                case CTSAR_DEEP:
                case CTSAR_START:
                    it = s->pass_iterations / s->depth;
                    if(s->iteration == s->depth - 1){
                        it = get_remaining(s);
                    }
                    if(s->iteration == 1 && s->type == CTSAR_START) {
                        it = s->pass_iterations - s->pass_iterations / s->depth;
                    }
                    break;
                default:
                    break;
            }
            calculate_multi_ratio(s,it,1);

            if(s->type == CTSAR_START && s->iteration > 0){
                s->type=CTSAR_DYNAMIC;
            }

            s->iteration++;//record loop iteration
        }
        pthread_barrier_wait(&s->all_barr);
        // s->all_barr->wait();
    }else{
        int dev=get_tid(),
            iters=s->chunk;

        if(s->old_devices[get_tid()].active){
            ctsar_device *d = &s->old_devices[get_tid()];
            ctsar_device *cd = &s->old_devices[s->num_devices-1];
            //ctsar_dev_group *g = &s->dev_groups[d->group];
            //ctsar_dev_group *cg = &s->dev_groups[cd->group];
            /* printf("chunk:%d, get_tid():%d\n", chunk, get_tid()); */
            if(ctsar_get_type(s)!=CTSAR_DEV_CPU){//decide how much to give the GPU
                if(s->type == CTSAR_CHUNK_STATIC){
                    iters = (1/s->ratio) * s->chunk;
                }else if (s->type == CTSAR_CHUNK_DYNAMIC | s->type == CTSAR_HYBRID_CHUNK){
                    // double cpi = s->dev_groups[s->old_devices[s->num_devices-1].group].avg_tpi;
                    if(d->pass_time == 0 || cd->pass_time == 0 || d->avg_tpi == 0){
                        iters = (1/s->ratio) * s->chunk;
                        /* iters = s->chunk; */
                    }else{
                        if(d->avg_tpi > d->tpi){//most recent was faster than average
                            if(d->requested_iterations <= d->iterations){//the change was due to having more iterations
                                /* if(get_tid()==0)
                                fprintf(stderr,"%d: increasing chunk size\n", get_tid()); */
                                iters = d->requested_iterations + s->chunk;//give it even more
                            }else{
                                iters = d->requested_iterations - s->chunk;//take even more back
                                /* if(get_tid()==0)
                                fprintf(stderr,"%d: decreasing chunk size\n", get_tid()); */
                            }
                        }else{//most recent was slower
                            if(d->requested_iterations <= d->iterations){//the change was due to having more iterations
                                iters = d->requested_iterations - s->chunk;//take even more back
                                /* if(get_tid()==0)
                                fprintf(stderr,"%d: decreasing chunk size\n", get_tid()); */
                            }else{
                                /* if(get_tid()==0)
                                fprintf(stderr,"%d: increasing chunk size\n", get_tid()); */
                                iters = d->requested_iterations + s->chunk;//give it even more
                            }
                        }
                    }

                    //HYBRID CHUNK is reset to adaptive in _loop
                    // DBG(DBG_TEST,"group[%d]: ratio=%lf, its=%d\n", get_tid(), rat, iters);
                }else{
                    iters = s->chunk;
                }
            }
            if(iters < 0){
                iters = 0;
            }
            s->old_devices[get_tid()].requested_iterations = iters;

            pthread_mutex_t assignment_mutex = PTHREAD_MUTEX_INITIALIZER;
            {
                pthread_mutex_lock(&assignment_mutex);
                assign_iterations(s,1, &dev, &iters);
                pthread_mutex_unlock(&assignment_mutex);
            }
        }
    }
    if(nested)
        ctsar_cpu_barrier(s);

    s->old_devices[get_tid()].total_inner++;
    s->old_devices[get_tid()].inner++;
    if(get_tid() == 0)
    {
        gettimeofday(&ctsar_time_end, NULL);
        uint64_t snt = calc_time(ctsar_time_start, ctsar_time_end);
        s->ctsar_next_time += snt;
        // DBG(DBG_TEST, "ctsar_next time, this pass=%lu, total=%lu\n", snt, s->ctsar_next_time);
    }
}

void  * ctsar_reg_reduc(ctsar * c, void * ptr, void * identity, size_t item_size, reduc_function reduc){
    int i;
    //void * buffer = NULL;
    ctsar_reduction * red = NULL;
    //ctsar_device *dev = &c->old_devices[get_tid()];
    for(i=0; i < c->num_reductions; i++){
        if(c->reductions[i].reduc_var == ptr){
            red = &c->reductions[i];
            break;
        }
    }
    if(red == NULL){
        //static ctsar_reduction *pass; //allows thread 0 to pass new reduction to other threads
        // c->all_barr->wait();
        pthread_barrier_wait(&c->all_barr);
        if(get_tid()==0)
        {
            DBG(DBG_DEBUG, "[%d] registering a new reduction region ptr=%p\n", get_tid(),ptr);
            if( red == NULL ){
                red = &c->reductions[c->num_reductions++];
            }
            red->item_size = item_size;
            red->reduc = reduc;
            red->reduc_var = ptr;
            red->identity = identity;

            red->host_data = ctsar_calloc_test(item_size, c->num_devices);
            /* fprintf(stderr,"host=%p\n", red->host_data); */

            char * source = (char*) identity;
            for(i=0; i<c->num_devices; i++){
                char * host_buff = (char*) red->host_data+i;
                for(int j=0; j<item_size; j++){
                    host_buff[j] = source[j];
                }
            }
        }
        // c->all_barr->wait();
        pthread_barrier_wait(&c->all_barr);
        for(i=0; i < c->num_reductions; i++){
            if(c->reductions[i].reduc_var == ptr){
                red = &c->reductions[i];
                break;
            }
        }
        if(red == NULL){
            fprintf(stderr, "this should be impossible, failing badly\n");
            exit(1);
        }
    }
    /* fprintf(stderr, "%d: registering %p offset %d from %p\n",get_tid(),((char*)red->host_data) + (get_tid() * item_size), get_tid()*item_size, red->host_data); */

    return ctsar_reg_mem(c, ((char*)red->host_data) + (get_tid() * item_size), item_size, 1,
            (ctsar_mem)(CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT));//TODO: consider releasing these
}

void complete_reductions(ctsar * c){
    //void * buffer = NULL;
    ctsar_reduction * red = NULL;
    //ctsar_device *dev = &c->old_devices[get_tid()];
    if(c->num_reductions > 0){
        pgi_sync(PGI_DEV_NOT_HOST);
        if(get_tid() == 0){
            for(int i=0; i < c->num_reductions; i++){
                red = &c->reductions[i];
                /* fprintf(stderr, "PERFORMING REDUCTION!\n"); */

                for(int j=0; j < c->num_devices; j++){
                    /* fprintf(stderr,"%d = %d\n", j, ((int*)red->host_data)[j]); */
                    red->reduc(red->reduc_var, red->reduc_var, ((char*)red->host_data) + (j * red->item_size));
                }
            }
        }
        // c->all_barr->wait();
        pthread_barrier_wait(&c->all_barr);
    }
}

int ctsar_loop(ctsar * s)
{
    int assigned;
    /* exit(1); */
    if(s->old_devices[get_tid()].active || (nested && allowed_devs & CTSAR_DEV_CPU) || s->type < CTSAR_CHUNK )//omp_get_thread_num() >= s->num_dev_groups && allowed_devs & ctsar_DEV_CPU))//not nice, but functional
    {

        if(s->type >= CTSAR_CHUNK && nested && s->old_devices[get_tid()].type == CTSAR_DEV_CPU){//to solve deadlock caused by GPU thread updating the fully_assigned counter while one or more threads is outside
            static int con=0;
            if(s->old_devices[get_tid()].active){
                con = s->fully_assigned;
            }
            ctsar_cpu_barrier(s);
            assigned = con;
        }else{
            assigned = s->fully_assigned;
        }
        if(assigned){
            /* if(sdebug > 2) fprintf(stderr,"%d testing loop, finishing\n",omp_get_thread_num()); */
            ctsar_device * dev = s->old_devices+(get_tid());
            //handle reductions
            if(dev->type != CTSAR_DEV_CPU){
                pgi_sync(PGI_DEV_NOT_HOST);
                for(int i=0; i < dev->num_regions; i++){
                    ctsar_region * reg = &dev->regions[i];
                    if( ! (reg->flags & CTSAR_MEM_PERSIST) )
                        reg->present = 0;
                }
            }
            if(s->type == CTSAR_HYBRID_CHUNK){
                s->type = CTSAR_DYNAMIC;
            }

            if(s->type >= CTSAR_CHUNK
                    && stop_slow
                    && dev->type == CTSAR_DEV_GPU
                    /* && s->type < CTSAR_CHUNK_DYNAMIC) */
                ){//adding stop slow to chunk schedulers
                    ctsar_device *d = dev;
                    ctsar_device *cd = &s->old_devices[s->num_devices-1];
                    ctsar_dev_group *g = &s->dev_groups[dev->group];
                    ctsar_dev_group *cg = &s->dev_groups[cd->group];
                    if(dev->pass_iterations <= cd->pass_iterations)//direct inner compares don't work with ratio chunk schedulers
                        ++(g->fall_count);
                    else
                        g->fall_count=0;
                    if(g->fall_count >= 2)
                    {
                        DBG(DBG_ALWAYS,"STOP_SLOW: converting group=%d to CPU threads\n",d->group);
                        /* ctsar_device * dev;
                           for(int j=0; j<g->size; j++){
                           dev = &s->old_devices[g->members[j]];
                           dev->tpi = cg->avg_tpi;
                           dev->type = CTSAR_DEV_CPU;
                           } */
                        d->tpi = cg->avg_tpi;
                        d->type = CTSAR_DEV_CPU;

                        g->avg_tpi = cg->avg_tpi;
                        g->type = CTSAR_DEV_CPU;

                        int cnt=0;
                        for(int j=0; j<s->num_devices;j++){
                            if(s->old_devices[j].type == CTSAR_DEV_CPU){
                                cnt++;
                            }
                        }
                        if(cnt > 0){
                            // s->b_barr = make_shared<boost::barrier>(cnt);
                            pthread_barrier_init(&s->b_barr, NULL, cnt);
                        }
                        //pthread_barrier_init(&s->barr, NULL, cnt);
                    }
                }

            ctsar_print_messages(s);

            complete_reductions(s);

            return 0;
        }else{
            /* if(sdebug > 2) fprintf(stderr,"%d testing loop, continuing\n",omp_get_thread_num()); */
            return 1;
        }
        /* return !s->fully_assigned; */
    }else{
        DBG(DBG_ALWAYS,"%d bailing, unneeded\n",pgi_get_thread_num());
        return 0;
    }
}

int ctsar_nested(){
    return nested;
    /* return 0; */
}
void ctsar_retarget_mem(ctsar * c, void * oldp, void * newp){
    int i;
    ctsar_device * dev = &c->old_devices[get_tid()];
    for(i=0; i < dev->num_regions; i++){
        if(dev->regions[i].host_data == oldp){
            dev->regions[i].host_data = newp;
        }
    }
}
void alloc_device_mem(ctsar_device * dev, ctsar_region *lreg){
     if(lreg->acc_data == nullptr){
        lreg->acc_data = pgi_malloc(lreg->device_pitch * lreg->rows * lreg->item_size);
    }
}
void *  ctsar_reg_mem_2d(ctsar * c, void * ptr, size_t item_size, size_t rows, size_t cols,
                      ctsar_mem flags, size_t pad_x, size_t pad_y, size_t *pitch_output){
    if(pitch_output != nullptr){
        *pitch_output = cols;
    }
    if(ctsar_get_type(c) == CTSAR_DEV_CPU)
        return ptr;

    int i;
    ctsar_device * dev = &c->old_devices[get_tid()];
    ctsar_region * reg = NULL;

    for(i=0; i < dev->num_regions; i++){
        ctsar_region * lreg = &dev->regions[i];
        if(lreg->host_data == ptr){
            if(       item_size != lreg->item_size
                    ||cols      != lreg->cols
                    ||rows      != lreg->rows){
                fprintf(stderr, "registered memory size change unimplemented\n");
                exit(1);
            }

            lreg->flags = flags;//allow reassigning of flags to reuse registrations
            /* DBG(-1,"returning preregistered memory\n"); */
            alloc_device_mem(dev, lreg);
            reg = lreg;
            goto return_ptr;
        }else if(lreg->host_data==NULL){
            reg = lreg;
        }
    }

    DBG(DBG_DEBUG, "[%d] registering a new region ptr=%p\n", get_tid(),ptr);
    if( reg == NULL ){
        reg = &dev->regions[dev->num_regions++];
    }

    reg->flags = flags;
    char * env;
    if(((env = getenv("CTSAR_NO_MEM_OPT")) != NULL) && atoi(env)){
        if(flags & CTSAR_MEM_OUTPUT){
            reg->flags = (ctsar_mem) (reg->flags & (CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT | CTSAR_MEM_INPUT));
        }else{
            reg->flags = (ctsar_mem) (reg->flags & (CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT));
        }
    }

    reg->host_data = ptr;
    reg->item_size = item_size;
    reg->rows = rows;
    reg->pad_x = pad_x;
    reg->pad_y = pad_y;
    reg->cols = cols;
    reg->device_pitch = cols;

    /* creg(ptr,item_size*count); */
    alloc_device_mem(dev, reg);
return_ptr:
        return reg->acc_data;
}
void *  ctsar_reg_mem(ctsar * c, void * ptr, size_t item_size, size_t count,
                      ctsar_mem flags){
    if(flags & CTSAR_MEM_COLUMN){
        fprintf(stderr, "column-wise partial copies require use of ctsar_reg_mem_2d\n");
        exit(EINVAL);
    }
    return ctsar_reg_mem_2d(c,ptr,item_size,count,1,flags,0,0, nullptr);//reuse 2d for 1d
}
void    ctsar_unreg_mem(ctsar * c, void * ptr){
    int i;
    ctsar_device * dev = &c->old_devices[get_tid()];
    for(i=0; i < dev->num_regions; i++){
        if(dev->regions[i].host_data == ptr || dev->regions[i].acc_data == ptr){
            pgi_free(dev->regions[i].acc_data);
            dev->regions[i].host_data = NULL;
            dev->regions[i].flags = (ctsar_mem)0;
            /* fprintf(stderr, "freed memory at address %p\n", ptr); */
        }
    }
}

void ctsar_swap_mem(ctsar * c, void * h1, void * h2){
    int i;
    ctsar_device * dev = &c->old_devices[get_tid()];
    for(i=0; i < dev->num_regions; i++){
        if(dev->regions[i].host_data == h1){
            dev->regions[i].host_data = h2;
            /* dev->regions[i].flags = 0; */
        }else{
            if(dev->regions[i].host_data == h2){
                dev->regions[i].host_data = h1;
            }
        }
    }
}

uint64_t ctsar_get_start(ctsar *c, int ltid){
    return c->old_devices[ltid].start;
}

uint64_t ctsar_get_end(ctsar *c, int ltid){
    return c->old_devices[ltid].end;
}

void * ctsar_calloc_test(size_t size, size_t count){
    void * ptr = NULL;

    ptr = calloc(size,count);

    if(ptr == NULL){
        perror("calloc failed");
        exit(errno);
    }
    return ptr;
}

void **ctsar_alloc_2d(size_t dts, size_t sa, size_t sb){
    void ** ret =  (void**)ctsar_calloc_test(sizeof(void *),sa);
    void * ptr = (void*)ctsar_calloc_test(dts,sa*sb);
    for(int i=0; i < sa; i++){
        ret[i] = (void*)(((char*)ptr) + ((i * sb) * dts));
    }
    return ret;
}

void ***ctsar_alloc_3d(size_t dts, size_t sa, size_t sb, size_t sc){
    void *** ret =  (void***)ctsar_calloc_test(sizeof(void **),sa);
    void * ptr = (void*)ctsar_calloc_test(dts,sa*sb*sc);
    for(int i=0; i < sa; ++i){
        ret[i] = (void **)ctsar_calloc_test(sizeof(void *), sb);
        for(int j=0; j < sb; ++j){
            ret[i][j] = (void *)((char *)ptr + (((sb * sc) * i) + (sc * j)) * dts);
        }
    }
    return ret;
}

};

