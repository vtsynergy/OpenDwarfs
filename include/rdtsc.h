#ifndef __RDTSC_H__
#define __RDTSC_H__
#include <time.h>
#define CHECK_ERROR(err) {if (err != CL_SUCCESS) { \
	fprintf(stderr, "Error1: %d\n", err);\
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


cl_event myEvent;
cl_ulong startTime, endTime;
unsigned long long int h2d_t = 0, k_t = 0, d2h_t = 0;
//bit field for faster ops
#ifdef ENABLE_TIMER
//#define TIMER_WORLD 0
#define TIMER_NAME_LEN 16

#define TIMER_D2H 1
#define TIMER_H2D 2
#define TIMER_D2D 4
#define TIMER_KERNEL 8
#define TIMER_USER 16
#define TIMER_COMPOSED -1 //use negative values for caomposed timers, so we can look at MSB as a quick identifier

enum timer_types {
    d2h = TIMER_D2H, h2d = TIMER_H2D, d2d = TIMER_D2D, kern = TIMER_KERNEL, user = TIMER_USER, comp = TIMER_COMPOSED
};

struct internalSingleTimer {
    enum timer_types type;
    char name[TIMER_NAME_LEN]; //Optional timer name, can be searched for grouping
    cl_ulong starttime, endtime;
    cl_event event;
};
//use the same struct except have an extra event

struct internalComposedTimer {
    enum timer_types type;
    char name[TIMER_NAME_LEN];
    cl_ulong starttime, endtime;
    cl_event event[2];
};

union internalTimer {
    struct internalSingleTimer s;
    struct internalComposedTimer c;
};

struct timer_group_mem;

struct timer_group_mem {
    union internalTimer * timer;
    struct timer_group_mem * next;
} head = {NULL, NULL}; //sentinel

struct timer_group_mem * tail;

//linear search for the timer

void * searchAll(union internalTimer * t) {
    //always skip sentinel
    struct timer_group_mem * curr = head.next;
    while (curr != 0 && curr->timer != t) {
        curr = curr->next;
    }
    if (curr != 0) {
        return (void *) curr;
    }
    return (void *) - 1;
}

int checkName(char * n1, char * n2, int len) {
    int i;
    if (n1 == 0 || n2 == 0)
        return -1;
    for (i = 0; i < len; i++) {
        if (n1[i] != n2[i]) return -1;
    }
    return 0;
}

void * searchAllByName(char * name, int len) {
    //always skip sentinel
    struct timer_group_mem * curr = head.next;
    while (curr != 0 && checkName(curr->timer->s.name, name, TIMER_NAME_LEN) != 0) {
        curr = curr->next;
    }
    if (curr != 0) {
        return (void *) curr;
    }
    return (void *) - 1;
}

//returns total time taken by all timers of a given type

cl_ulong aggregateTimesFromType(int type) {
    //always skip sentinel
    struct timer_group_mem * curr = head.next;
    cl_ulong time = 0;
    while (curr != 0 && curr->timer != 0) {
        if (curr->timer->s.type == type && curr->timer->s.endtime >= curr->timer->s.starttime) {
            time += curr->timer->s.endtime - curr->timer->s.starttime;
        }
        curr = curr->next;
    }
    return time;
}
//Adds a new timer_group object to the groups list, and returns its integer ID

//int newTimerGroup() {
//timer_group curr = timer_world;
//	while (curr.next != 0) {
//		curr = *curr.next;
//	}
//curr.next = (struct timer_group *) malloc(sizeof(struct timer_group), 1);
//curr.next->groupID = curr.groupID+1;
//return curr.groupID+1;
//}

//only returns the primary timer, not any composed timers
void * getTimePtr(cl_event e) {
    struct timer_group_mem * curr = head.next;
    while (curr != 0) {
        //if composed, will be a negative value
        if (curr->timer->s.type > 0) {
            if (curr->timer->s.event == e) return (void *) curr->timer;
        }
        curr = curr->next;
    }
    return (void *) -1;
}
//only returns a composed timer with events matching both e1 and e2, in either order
void * getDualTimePtr(cl_event e1, cl_event e2) {
    struct timer_group_mem * curr = head.next;
    while (curr != 0) {
        //if composed, will be a negative value
        if(curr->timer->s.type < 0) {
            if((curr->timer->c.event[0] == e1 \
                    && curr->timer->c.event[1] == e2) \
                    || (curr->timer->c.event[0] == e2 \
                    && curr->timer->c.event[1] == e1)) \
                return (void *) curr->timer;
        }
        curr = curr->next;
    }
    return (void *) -1;
}
void walkList() {
    struct timer_group_mem * curr = head.next;
    printf("Walking list starting at [%lx]--[%lx]\n", (unsigned long) &head, (unsigned long) head.timer);
    while (curr != 0) {
        printf("\t[%lx]--[%lx]\n", (unsigned long) curr, (unsigned long) curr->timer->s.event);
        curr = curr->next;
    }
}
//Adds the timer to the specified group
//Updates group structure and timer
//int addToTimerGroup(struct internalTimer * timer, timer_group *group) {
//if (group->groupID < 6) { //Don't let the user manually edit core groups
//fprintf(stderr, "Error: Groups 0-5 cannot be manually edited.\n");
//return -1;
//}

//searchGroup(g, t)
//}

//Removes the timer from the specified group
//Updates group structure and timer
//int removeFromTimerGroup(struct internalTimer * timer, timer_group *group) {
//if (group->groupID < 6) { //Don't let the user manually edit core groups
//fprintf(stderr, "Error: Groups 0-5 cannot be manually edited.\n");
//return -1;
//}

//}

//simply adds timer t to the end of the list

addTimer(union internalTimer * t) {
    if (head.next == NULL) { //no members
        tail = &head; //reset tail, just incase
    } else {
        while (tail->next != NULL) { //slide tail to the end
            tail = tail->next;
        }
    }
    struct timer_group_mem * temp_wrap = (struct timer_group_mem *) malloc(sizeof (struct timer_group_mem));

    temp_wrap->next = NULL;
    temp_wrap->timer = t;
    tail->next = temp_wrap;
    tail = temp_wrap;
}

//irreversible! Only call immediately before freeing the timer!

int removeTimer(union internalTimer * t) {
    struct timer_group_mem * curr = head.next, * old = &head;
    while (curr != 0 && curr->timer != t) {
        old = curr;
        curr = curr->next;
    }
    if (curr != 0) {
        if (curr->next == 0) { //we are the tail!
            tail = old; //so back the tail up one
        }
        old->next = curr->next;
        free(curr);
        return 0;
    }
    return -1; // probably should free(groups) somehow 
}

//should work for composed timers, so long as start_timer is strictly used on the first of the two events
#define START_TIMER(e) {void * ptr = getTimePtr(e); \
                        if (ptr == (void *) -1) {\
                                fprintf(stderr, "Timer Error: Cannot start uninitialized timer for event [%lx]!\n", (unsigned long) e); \
                        } else {\
                            cl_int err = clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &((union internalTimer *)ptr)->s.starttime, NULL); \
                            CHECK_ERROR(err)\
                        }\
        }
#define INI_TIMER(e, t) {void * ptr = getTimePtr(e); \
                         if (ptr == (void *) -1) {\
                                 if(t >= TIMER_USER || t <= TIMER_COMPOSED) { \
                                         fprintf(stderr, "Timer Error: invalid type [%d] for INI_TIMER!\nTimer for event [%lx] not initialized!", t, (unsigned long) e); }\
                                 struct internalSingleTimer * temp = (struct internalSingleTimer*) malloc(sizeof(struct internalSingleTimer)); \
                                 temp->type = t;\
                                 temp->starttime = 0;\
                                 temp->endtime = 0;\
                                 temp->event = e;\
                                 addTimer((union internalTimer *)temp);\
                         } else { \
                 		fprintf(stderr, "Timer Error: event [%lx] has already initialized timer [%lx]!\n", (unsigned long) e, (unsigned long) ptr); \
                  	}}
#define DEST_TIMER(e) {void * ptr = getTimePtr(e); \
                       if(ptr != (void *) -1) { \
                       		removeTimer((union internalTimer *) ptr); \
                       		free(ptr); \
                       }}
//should work for composed timers, so long as end_timer is strictly used on the last
#define END_TIMER(e) {void * ptr = getTimePtr(e); \
                        if (ptr == (void *) -1) {\
                                fprintf(stderr, "Timer Error: Cannot end uninitialized timer for event [%lx]!\n", (unsigned long) e); \
                        } else { \
                            cl_int err = clGetEventProfilingInfo(e, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &((union internalTimer *)ptr)->s.endtime, NULL); \
                            CHECK_ERROR(err)\
                        }}
#define TIMER_ENABLE CL_QUEUE_PROFILING_ENABLE
#define TOTAL_H2D aggregateTimesFromType(TIMER_H2D)
#define TOTAL_D2H aggregateTimesFromType(TIMER_D2H)
#define TOTAL_D2D aggregateTimesFromType(TIMER_D2D)
#define TOTAL_KERNEL aggregateTimesFromType(TIMER_KERNEL)
#define TOTAL_USER aggregateTimesFromType(TIMER_USER)
#define TOTAL_COMPOSED aggregateTimesFromType(TIMER_COMPOSED)
#define PRINT_CORE_TIMERS {printf("********************************************************************************\n"\
"OCD Core Timers\n"\
"********************************************************************************\n"\
"\tHost to Device Time:   [%lu]\n"\
"\tDevice to Host Time:   [%lu]\n"\
"\tDevice to Device Time: [%lu]\n"\
"\tDevice Kernel Time:    [%lu]\n"\
"\tUser Timer Total:      [%lu]\n"\
"\tComposed Timer Total:  [%lu]\n"\
"********************************************************************************\n"\
,TOTAL_H2D, TOTAL_D2H, TOTAL_D2D, TOTAL_KERNEL, TOTAL_USER, TOTAL_COMPOSED);}
#define TIMER_DEBUG_WALK_LIST {walkList();}

//#define START_DUAL_TIMER(a, b)
//#define INI_DUAL_TIMER(a, b)
//#define DEST_DUAL_TIMER(a, b)
//#define END_DUAL_TIMER(a, b)
//#define PRINT_COUNT printf("H2D Transfer Time: %f (ms)\n",h2d_t*1e-6f);\
		    printf("Kernel Time: %f (ms)\n",k_t*1e-6f);\
		    printf("D2H Transfer Time: %f (ms)\n",d2h_t*1e-6f);

#define CL_FINISH(c) clFinish(c);
#else
//#define INI_TIMER(e, t)
#define INI_TIMER
//#define DEST_TIMER(e)
//#define START_TIMER(e)
#define START_TIMER
#define TIMER_ENABLE 0
//#define END_TIMER(e)
#define END_TIMER
#define COUNT_H2D
#define COUNT_K
#define COUNT_D2H
#define PRINT_COUNT 
#define CL_FINISH(c)
#endif

#ifdef START_POWER
#define START_KERNEL printf("Kernel Start\n");
#define END_KERNEL printf("Kernel END\n");
#else
#define START_KERNEL
#define END_KERNEL
#endif

cl_device_id
GetDevice(int platform, int device) {
    cl_int err;
    cl_uint nPlatforms = 1;
    err = clGetPlatformIDs(0, NULL, &nPlatforms);
    CHECK_ERROR(err);

    if (nPlatforms <= 0) {
        printf("No OpenCL platforms found. Exiting.\n");
        exit(0);
    }
    if (platform < 0 || platform >= nPlatforms) // platform ID out of range
    {
        printf("Platform index %d is out of range. \n", platform);
        exit(-4);
    }
    cl_platform_id *platforms = (cl_platform_id *) malloc(sizeof (cl_platform_id) * nPlatforms);
    err = clGetPlatformIDs(nPlatforms, platforms, NULL);
    CHECK_ERROR(err);

    cl_uint nDevices = 1;
    char platformName[100];
    err = clGetPlatformInfo(platforms[0], CL_PLATFORM_VENDOR, sizeof (platformName), platformName, NULL);
    CHECK_ERROR(err);
    printf("Platform Chosen : %s\n", platformName);
    // query devices
    err = clGetDeviceIDs(platforms[platform], USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 0, NULL, &nDevices);
    CHECK_ERROR(err);
    if (nDevices <= 0) {
        printf("No OpenCL Device found. Exiting.\n");
        exit(0);
    }
    if (device < 0 || device >= nDevices) // platform ID out of range
    {
        printf("Device index %d is out of range. \n", device);
        exit(-4);
    }
    cl_device_id* devices = (cl_device_id *) malloc(sizeof (cl_device_id) * nDevices);
    err = clGetDeviceIDs(platforms[platform], USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, nDevices, devices, NULL);
    char DeviceName[100];
    err = clGetDeviceInfo(devices[device], CL_DEVICE_NAME, sizeof (DeviceName), DeviceName, NULL);
    CHECK_ERROR(err);
    printf("Device Chosen : %s\n", DeviceName);

    return devices[device];
}


#endif

