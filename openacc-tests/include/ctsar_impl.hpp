#ifndef __CTSAR_IMPL
#define __CTSAR_IMPL

#include "config.h"
#include <ctsar.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <memory>
#include <pthread.h>
// #include <boost/thread/barrier.hpp>
// #include <boost/format.hpp>
#include <hwloc.h>
// #include "matrix.hpp"


class ctsar{
  public:
    ctsar();
    ~ctsar();

    ctsar_type type;
    pthread_barrier_t all_barr;
    pthread_barrier_t b_barr;

    int num_devices;//no longer specifically treating as CPU or GPU
    ctsar_device * old_devices;//indexed by thread number in independent mode, split in nested mode
    // std::vector<std::shared_ptr<CTSAR::Device>> devices_;

    int num_dev_groups;//the number of schedulable entities for use in the LP model
    ctsar_dev_group * dev_groups;

    int fully_assigned;

    size_t pass_iterations;
    size_t assigned;  //number of iterations assigned so far

    int depth;
    int chunk;
    double ratio;

    int fallback; //boolean, whether to fall back on basic omp dynamic schedule

    size_t iteration;
    size_t outer_iteration;

    //global aggregators for messages
    uint64_t ctsar_next_time;

    ctsar_reduction reductions[5]; //TODO: unsafe, but expedient
    size_t num_reductions;
};

#endif //__CTSAR_IMPL
