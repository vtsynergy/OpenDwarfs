#ifndef __RDTSC_H__
#define __RDTSC_H__

#ifdef __cplusplus
extern "C" {
#endif


#include <sys/time.h>
#define CHECK_ERROR(err) {if (err != CL_SUCCESS) { \
	fprintf(stderr, "Error: %d\n", err);\
	exit(1); \
}}
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#define USEGPU 1
#define PLATFORM_ID 0
#define DEVICE_ID 0
#include <stdint.h>
#ifndef _STRING_H
#include <string.h>
#endif
#include <stdio.h>


#include "../build/config.h"



extern cl_event ocdTempEvent;

#ifdef ENABLE_TIMER
//use negative values for composed timers, so we can potentially look at MSB as a quick identifier
enum timer_types {
    OCD_TIMER_D2H = 1, OCD_TIMER_H2D = 2, OCD_TIMER_D2D = 4, OCD_TIMER_KERNEL = 8, OCD_TIMER_HOST = 16, OCD_TIMER_DUAL = -1
};

struct ocdTimer {
    enum timer_types type;
    const char * name;
    int nlen;
    //char name[TIMER_NAME_LEN]; //Optional timer name, can be searched for grouping
    cl_ulong starttime, endtime;
    cl_event event;
};
//use the same struct except have an extra event


extern struct ocdTimer * ocdTempTimer;

struct ocdDualTimer {
    enum timer_types type;
    const char * name;
    int nlen;
    //char name[TIMER_NAME_LEN];
    cl_ulong starttime, endtime;
    cl_event event[2];
};


extern struct ocdDualTimer * ocdTempDualTimer;


//host timers don't actually use events, rather two gettimeofday calls
// which return time values immediately
//microsecond resolution is scaled up by 1000 to be compatible with
//CL-based timers
struct ocdHostTimer {
    enum timer_types type;
    const char * name;
    int nlen;
    cl_ulong starttime, endtime;
    struct timeval timer;
} ;

extern struct ocdHostTimer fullExecTimer;
//the above fullExecTimer is a special purpose timer which does not reside on
//any list, but measures host time from OCD_INIT to OCD_FINISH

extern struct ocdHostTimer * ocdTempHostTimer;

#ifdef TIMER_TEST
//number of fake names to generate
#define TIMER_TEST_NAME_COUNT 200
//maximum number of fake events per name
#define TIMER_TEST_MAX_LIFE 20
#define TIMER_TEST_MAX_LENG 30
#endif



union ocdInternalTimer {
    struct ocdTimer s;
    struct ocdDualTimer c;
    struct ocdHostTimer h;
};

struct timer_group_mem;

struct timer_group_mem {
    union ocdInternalTimer * timer;
    struct timer_group_mem * next;
    struct timer_group_mem * alphanext; //ignored except for alpha sort
};
extern struct timer_group_mem head; //sentinel

extern struct timer_group_mem * tail;

extern char rootStr[1];
extern cl_ulong rootTimes[7]; 
extern cl_ulong totalTimes[7];

struct timer_name_tree_node {
    const char * string; //the first character is hijacked as a flag for pointer ownership
    //to make sense of who is responsible for freeing at the end
    //this lets one descendant branch reuse our space
    int len; //length, not counting flag and zero-byte
    struct timer_name_tree_node * next;
    struct timer_name_tree_node * child;
    struct timer_group_mem * n_head; //first list node for a timer matching this name
    int tcount;
    cl_ulong * times; //pointer to a 7-member array of cl_ulongs
    //one aggregator for each type, and another for all
}; 

extern struct timer_name_tree_node root;






//linear search of the Name List.
//returns a pointer to the correct time array, or -1 if none exists yet
//rather inefficient if many names are used, but the tree will take care of
// speeding lookups, and we'll switch to alpha sort by default as a sideffect
extern void * checkSimpleNameList(const char * s, int len);

extern struct timer_name_tree_node * atail;

//simple named timer aggregation
//linear scan of the timer list, adds nodes to a names list as necessary
//DO NOT USE AT THE SAME TIME AS THE TREE
//this replaces the tree with a simple unordered list
extern void simpleNameTally(); 

//assumes simpleNameTally was already called (once) to add up timers
//now culls off zero-value timers
extern void simpleNamePrint();

//chews up the timer list from head to tail, deallocating all nodes
extern void destTimerList();

//chews up the simpleNameList from root to atail, deallocating all nodes
extern void destNameList();

//only returns the primary timer, not any composed timers
extern void * getTimePtr(cl_event e); 

//only returns a composed timer with events matching both e1 and e2, in either order
extern void * getDualTimePtr(cl_event e1, cl_event e2);

//simply adds timer t to the end of the list
extern void addTimer(union ocdInternalTimer * t);

//irreversible! Only call immediately before freeing the timer!
extern int removeTimer(union ocdInternalTimer * t);


#ifdef TIMER_TEST
//Debug call for checking list construction
extern void walkList();

#endif



//should work for composed timers, so long as start_timer is strictly used on the first of the two events
#define START_TIMER(e, t, n, p) {void * ptr = getTimePtr(e); \
        if (ptr == (void *) -1) {\
                /*fprintf(stderr, "Timer Error: Cannot start uninitialized timer for event [%lx]!\n", (unsigned long) e); */\
                if(t >= OCD_TIMER_HOST || t <= OCD_TIMER_DUAL) { \
                        fprintf(stderr, "Timer Error: invalid type [%d] for START_TIMER!\nTimer for event [%lx] not initialized or started!", t, (unsigned long) e); \
                }else { \
                        struct ocdTimer * temp = (struct ocdTimer*) calloc(sizeof(struct ocdTimer), 1); \
                        temp->type = t;\
                        temp->event = e;\
                            temp->name = n;\
                        addTimer((union ocdInternalTimer *)temp);\
                        p = temp; /*set the return pointer*/\
                        cl_int err = clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof (cl_ulong), &temp->starttime, NULL); \
                        CHECK_ERROR(err)\
                }\
        } else {\
                p = &((union ocdInternalTimer *)ptr)->s; /*set the return pointer*/\
                cl_int err= clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof (cl_ulong), &((union ocdInternalTimer *)ptr)->s.starttime, NULL); \
                CHECK_ERROR(err)\
        }\
}

//starts a gettimeofday-based timer
#define START_HOST_TIMER(n, p) {\
struct ocdHostTimer * temp = (struct ocdHostTimer*) calloc(sizeof(struct ocdHostTimer), 1);\
temp->type = OCD_TIMER_HOST;\
gettimeofday(temp->timer, NULL);\
temp->starttime = 1000 * (temp->timer.tv_sec*1000000L + temp->timer.tv_usec);\
p = temp;\
}

//should work for composed timers, so long as end_timer is strictly used on the last
//assumes t is a valid single event, does not check if it's a dual
#define END_TIMER(t) {\
                            cl_int err = clGetEventProfilingInfo(t->event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), (void *)&t->endtime, NULL); \
                            CHECK_ERROR(err)\
                        }

//assumes t is a valid timer, ensures it's a host-type
#define END_HOST_TIMER(t) {\
if (t->type == OCD_TIMER_HOST) {\
gettimeofday(t->timer, NULL);\
temp->endtime = 1000 * (t->timer.tv_sec*1000000L + t->timer.tv_usec);\
}\
}
#define TOTAL_EXEC totalTimes[0]
#define TOTAL_D2H totalTimes[1]
#define TOTAL_H2D totalTimes[2]
#define TOTAL_D2D totalTimes[3]
#define TOTAL_KERNEL totalTimes[4]
#define TOTAL_HOST totalTimes[5]
#define TOTAL_DUAL totalTimes[6]
#define OCD_PRINT_TIMERS {printf("********************************************************************************\n"\
"OCD Core Timers (nanoseconds)\n"\
"********************************************************************************\n"\
"Total Execution Time:  \t[%lu]\n"\
"\tHost to Device Time:   [%lu]\n"\
"\tDevice to Host Time:   [%lu]\n"\
"\tDevice to Device Time: [%lu]\n"\
"\tDevice Kernel Time:    [%lu]\n"\
"\tUser Timer Total:      [%lu]\n"\
"\tComposed Timer Total:  [%lu]\n"\
"********************************************************************************\n"\
,TOTAL_EXEC, TOTAL_H2D, TOTAL_D2H, TOTAL_D2D, TOTAL_KERNEL, TOTAL_HOST, TOTAL_DUAL);}


//absolutely everything needed to start the timers
#define TIMER_INIT {\
gettimeofday(&fullExecTimer.timer, NULL);\
fullExecTimer.starttime = 1000 * (fullExecTimer.timer.tv_sec*1000000L + fullExecTimer.timer.tv_usec);\
}

//and absolutely everything needed to finalize them
// performs timer aggregation and printing
// deconstructs timer list and name tree/list
    //TODO-free all our data structures
#define TIMER_FINISH {\
gettimeofday(&fullExecTimer.timer, NULL);\
fullExecTimer.endtime = 1000 * (fullExecTimer.timer.tv_sec*1000000L + fullExecTimer.timer.tv_usec);\
    simpleNameTally();\
    OCD_PRINT_TIMERS\
    simpleNamePrint();\
    destNameList();\
    destTimerList();\
    }
//starts the dual timer specified by events a and b, assumes a is the "first" event
#define START_DUAL_TIMER(a, b, n, p) {void * ptr = getDualTimePtr(a, b); \
                        if (ptr == (void *) -1) {\
                                /*fprintf(stderr, "Timer Error: Cannot start uninitialized timer for events [%lx] and [%lx]!\n", (unsigned long) a, (unsigned long) b);*/ \
                                 struct ocdDualTimer * temp = (struct ocdDualTimer*) calloc(sizeof(struct ocdDualTimer), 1); \
                                 temp->type = OCD_TIMER_DUAL;\
                                 temp->event[0] = a;\
                                 temp->event[1] = b;\
                            temp->name = n;\
                                 addTimer((union ocdInternalTimer *)temp);\
                                 p = temp; /*set the return pointer*/\
                                 cl_int err = clGetEventProfilingInfo(a, CL_PROFILING_COMMAND_START, sizeof (cl_ulong), &temp->starttime, NULL); \
                                 CHECK_ERROR(err)\
                        } else {\
                            p = &((union ocdInternalTimer *)ptr)->c;\
    cl_int err = clGetEventProfilingInfo(a, CL_PROFILING_COMMAND_START, sizeof (cl_ulong), &((union ocdInternalTimer *)ptr)->s.starttime, NULL); \
                            CHECK_ERROR(err)\
                        }\
        }
//assumes t is a valid ocdTimer, but ensures it's a dual timer
#define END_DUAL_TIMER(t) {\
                        if (t->type == OCD_TIMER_DUAL) { \
                            cl_int err = clGetEventProfilingInfo(t->event[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), (void *)&t->endtime, NULL); \
                            CHECK_ERROR(err)\
                        }}

#else
#define OCD_INIT
#define OCD_FINISH
#define START_TIMER(e, t, n, p)
#define END_TIMER(t)
#define START_DUAL_TIMER(a, b, n, p)
#define END_DUAL_TIMER(t)
#define PRINT_CORE_TIMERS
#define START_HOST_TIMER(n, p)
#define END_HOST_TIMER(t)
#endif

#ifdef START_POWER
#define START_KERNEL printf("Kernel Start\n");
#define END_KERNEL printf("Kernel END\n");
#else
#define START_KERNEL
#define END_KERNEL
#endif

extern cl_device_id GetDevice(int platform, int device);

#ifdef __cplusplus
}
#endif

#endif //FILE
