#ifndef __RDTSC_H__
#define __RDTSC_H__
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



//#define ENABLE_TIMER


cl_event ocdTempEvent;
cl_ulong startTime, endTime;
unsigned long long int h2d_t = 0, k_t = 0, d2h_t = 0;
//use negative values for composed timers, so we can look at MSB as a quick identifier
enum timer_types {
    OCD_TIMER_D2H = 1, OCD_TIMER_H2D = 2, OCD_TIMER_D2D = 4, OCD_TIMER_KERNEL = 8, OCD_TIMER_HOST = 16, OCD_TIMER_DUAL = -1
};

struct ocdTimer {
    enum timer_types type;
    char * name;
    int nlen;
    //char name[TIMER_NAME_LEN]; //Optional timer name, can be searched for grouping
    cl_ulong starttime, endtime;
    cl_event event;
};
//use the same struct except have an extra event


struct ocdTimer * ocdTempTimer;

struct ocdDualTimer {
    enum timer_types type;
    char * name;
    int nlen;
    //char name[TIMER_NAME_LEN];
    cl_ulong starttime, endtime;
    cl_event event[2];
};


struct ocdDualTimer * ocdTempDualTimer;


//host timers don't actually use events, rather two gettimeofday calls
// which return time values immediately
//microsecond resolution is scaled up by 1000 to be compatible with
//CL-based timers
struct ocdHostTimer {
    enum timer_types type;
    char * name;
    int nlen;
    cl_ulong starttime, endtime;
    struct timeval timer;
} fullExecTimer = {OCD_TIMER_HOST, NULL, 0, 0, 0, {0, 0}};
//the above fullExecTimer is a special purpose timer which does not reside on
//any list, but measures host time from OCD_INIT to OCD_FINISH

struct ocdHostTimer * ocdTempHostTimer;
//bit field for faster ops
#ifdef ENABLE_TIMER

#define TIMER_TEST
#ifdef TIMER_TEST
//number of fake names to generate
#define TIMER_TEST_NAME_COUNT 200
//maximum number of fake events per name
#define TIMER_TEST_MAX_LIFE 20
#define TIMER_TEST_MAX_LENG 30
#endif
//#define TIMER_WORLD 0
#define TIMER_NAME_LEN 16



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
} head = {NULL, NULL, NULL}; //sentinel

struct timer_group_mem * tail;

char rootStr[2] = {(char)1, (char)0};
cl_ulong rootTimes[7] = {0, 0, 0, 0, 0, 0, 0};
cl_ulong totalTimes[7] = {0, 0, 0, 0, 0, 0, 0};

struct timer_name_tree_node {
    char * string; //the first character is hijacked as a flag for pointer ownership
    //to make sense of who is responsible for freeing at the end
    //this lets one descendant branch reuse our space
    int len; //length, not counting flag and zero-byte
    struct timer_name_tree_node * next;
    struct timer_name_tree_node * child;
    struct timer_group_mem * n_head; //first list node for a timer matching this name
    int tcount;
    cl_ulong * times; //pointer to a 7-member array of cl_ulongs
    //one aggregator for each type, and another for all
} root = {
    rootStr, 0, NULL, NULL, &head, 0, rootTimes
}; //sentinel

//linear search for the timer
//
//void * searchAll(union ocdInternalTimer * t) {
//    //always skip sentinel
//    struct timer_group_mem * curr = head.next;
//    while (curr != 0 && curr->timer != t) {
//        curr = curr->next;
//    }
//    if (curr != 0) {
//        return (void *) curr;
//    }
//    return (void *) - 1;
//}

//returns 0 for a perfect match
//-1 if either pointer is null
//or the positive int index of the mismatched character (including null terminator)

int checkName(char * n1, char * n2, int len) {
    int i;
    if (n1 == 0 || n2 == 0)
        return -1;
    for (i = 1; i <= len + 1; i++) {
        if (n1[i] != n2[i]) return i;
    }
    return 0;
}

//descends the name tree as far as possible, either until an exact match or
// nearest insertion point is found, and returns a pointer to it
/*
struct timer_name_tree_node * descendTree(char*str, struct timer_name_tree_node * r, int slen) {
    struct timer_name_tree_node * node = r;
    int idx = 0;
    str -= 1; //back up to fake the flag byte
    while (node != 0)
        int cmp = checkName(node->string, str, node->len);
    if (cmp == 0) { //exact match of slen characters
        if (node->child != NULL) {
            if (node->child->string[0] <= str[cmp]) { //They're before me so descend
                str += cmp;
                slen -= cmp;
                node = node->child;
            } else {

            }
        }
        //else perfect match
        return node;

    } else if (cmp >= 2) { //prefix match
        if (node->len == cmp - 1) { //the stored node is a prefix to me
            if (node->child == NULL) { //not yet branched
                return node;
            }
            //else descend

            //else need to change parent's child
            return node;
        } else { //must branch
            return node;
        }

    } else if (cmp == 1) { //one string ran out
        //look ahead
        if (node)
        }
}
}*/

//traverses the tree, adding node branches as necessary
//returns a pointer to the perfect match node (new or otherwise)
//CURRENTLY UNUSED, UNFINISHED, AND UNTESTED
//struct timer_name_tree_node * addName(char*str, struct timer_name_tree_node * r, int slen) {
//    struct timer_name_tree_node * node = r;
//    int idx = 0;
//    str -= 1; //back up to fake the flag byte
//    while (node != 0) {
//        int cmp = checkName(node->string, str, node->len);
//        if (cmp == 0) { //exact match of len characters
//            if (node->child != NULL) {
//                //update pointers/indices
//                str += node->len;
//                slen -= node->len;
//                //if (node->child->string[1] <= str[1]) { 
//                //since it matched the extra tail character, we must have once existed
//                //so it's safe to just descend
//                node = node->child;
//                //} else {
//
//                // }
//            } else { //this is a perfect match leaf since it matched the null-terminator
//                return node;
//                //add a timer to the tail pointed to in node
//
//            }
//        } else if (cmp > 0) { //prefix match
//            if (node->len == cmp - 1) { //the stored node is a prefix to me
//                //update pointers/indices
//                str += cmp;
//                slen -= cmp;
//                if (node->child == NULL) { //not yet branched
//                    //create a new node for the rest of my string
//                    struct timer_name_tree_node * n = (struct timer_name_tree_node *) malloc(sizeof (struct timer_name_tree_node));
//                    n->next = NULL;
//                    n->len = slen;
//                    n->child = NULL;
//                    n->string = (char *) calloc(sizeof (char), slen + 2);
//                    memcpy(&n->string[1], &str[1], slen);
//                    //point the current node to a new child sentinel
//                    node->child = (struct timer_name_tree_node *) malloc(sizeof (struct timer_name_tree_node));
//                    node->child->next = n;
//                    node->child->len = 0;
//                    node->child->child = NULL;
//                    node->child->string = &node->string[cmp];
//                    //return the new node
//                    return n;
//                } else {
//                    //else descend
//                    node = node->child;
//                }
//            } else if (slen == cmp - 1) { //i should be a prefix to the stored node
//                //do something like the above, but reverse the order
//                //if (node->child == NULL) { //not yet branched
//                //create a new node for the sentinel
//                struct timer_name_tree_node * n1 = (struct timer_name_tree_node *) malloc(sizeof (struct timer_name_tree_node));
//                n1->len = 0;
//                n1->child = NULL;
//                n1->string = (char *) calloc(sizeof (char), 2);
//                memcpy(&n1->string[1], &str[1], slen);
//                //point the current node to a new child for the rest of my string
//                struct timer_name_tree_node * n2 = (struct timer_name_tree_node *) malloc(sizeof (struct timer_name_tree_node));
//                n1->next = n2;
//                n2->next = NULL;
//                n2->len = node->len - cmp + 1;
//                n2->child = node->child; //don't care whether it has a child or not, branch is upstream
//                n2->string = &node->string[cmp];
//                node->child = n2;
//                //shrink the already existing node
//                node->len = cmp - 1;
//                //} else {
//                //else descend
//                //}
//            } else { //we have something in common, but we need to split in the middle
//                //if (node->child == NULL) { //not yet branched
//                //create a node for the new branch
//                struct timer_name_tree_node * n1 = (struct timer_name_tree_node *) malloc(sizeof (struct timer_name_tree_node));
//                n1->len = slen - cmp + 1;
//                n1->child = NULL;
//                n1->string = (char *) calloc(sizeof (char), 1 + slen - cmp);
//                memcpy(&n1->string[1], &str[1], slen);
//                //point the current node to a new child for the rest of my string
//                struct timer_name_tree_node * n2 = (struct timer_name_tree_node *) malloc(sizeof (struct timer_name_tree_node));
//                n2->len = node->len - cmp + 1;
//                n2->child = node->child;
//                n2->string = &node->string[cmp];
//                //shorten the parent
//                node->len = cmp - 1;
//                if (str[1] < node->string[1 + cmp]) { //enforce alphabetic node ordering
//                    n1->next = n2;
//                    n2->next = NULL;
//                    node->child = n1;
//                } else {
//                    n2->next = n1;
//                    n1->next = NULL;
//                    node->child = n2;
//                }
//                // } else {
//                //else descend
//                // }
//            }
//
//
//        }
//    }
//}

//linear search of the Name List.
//returns a pointer to the correct time array, or -1 if none exists yet
//rather inefficient if many names are used, but the tree will take care of
// speeding lookups, and we'll switch to alpha sort by default as a sideffect
void * checkSimpleNameList(char * s, int len) {
    struct timer_name_tree_node * curr = root.next;
    while (curr != NULL) { //still unique names to be checked
        if (checkName(s, curr->string, len) == 0) {
            return curr->times;
        }
        curr = curr->next;
    }
    return (void *)-1;
}

struct timer_name_tree_node * atail = &root;
//simple named timer aggregation
//linear scan of the timer list, adds nodes to a names list as necessary
//DO NOT USE AT THE SAME TIME AS THE TREE
//this replaces the tree with a simple unordered list
void simpleNameTally() {
    struct timer_group_mem * curr = head.next;
    while (curr != NULL) {
        void * time;
        if (curr->timer->s.name != NULL) {
            time = checkSimpleNameList(curr->timer->s.name, curr->timer->s.nlen);
            if (time == (void *)-1) {
                //initialize a new name list node
                atail->next = (struct timer_name_tree_node *) calloc(sizeof (struct timer_name_tree_node), 1);
                atail = atail->next;
                atail->next = NULL;
                atail->len = curr->timer->s.nlen;
                atail->string = curr->timer->s.name;
                atail->times = (cl_ulong *) calloc(sizeof(cl_ulong), 7);
                time = (void *)atail->times;
            }} else {
            time = (void *)root.times;
            }
            if (curr->timer->s.endtime > curr->timer->s.starttime) {
                switch (curr->timer->s.type) {
                    case OCD_TIMER_D2H:
                ((cl_ulong *) time)[1] += curr->timer->s.endtime - curr->timer->s.starttime;
                totalTimes[1] +=curr->timer->s.endtime - curr->timer->s.starttime;
                break;
                
                    case OCD_TIMER_H2D:
                ((cl_ulong *) time)[2] += curr->timer->s.endtime - curr->timer->s.starttime;
                totalTimes[2] +=curr->timer->s.endtime - curr->timer->s.starttime;
                break;
                                        
                    case OCD_TIMER_D2D:
                ((cl_ulong *) time)[3] += curr->timer->s.endtime - curr->timer->s.starttime;
                totalTimes[3] +=curr->timer->s.endtime - curr->timer->s.starttime;
                break;
                        
                    case OCD_TIMER_KERNEL:
                ((cl_ulong *) time)[4] += curr->timer->s.endtime - curr->timer->s.starttime;
                totalTimes[4] +=curr->timer->s.endtime - curr->timer->s.starttime;
                break;
                        
                    case OCD_TIMER_HOST:
                ((cl_ulong *) time)[5] += curr->timer->s.endtime - curr->timer->s.starttime;
                totalTimes[5] +=curr->timer->s.endtime - curr->timer->s.starttime;
                break;
                        
                    case OCD_TIMER_DUAL:
                ((cl_ulong *) time)[6] += curr->timer->s.endtime - curr->timer->s.starttime;
                totalTimes[6] +=curr->timer->s.endtime - curr->timer->s.starttime;
                break;
                }
                ((cl_ulong *) time)[0] += curr->timer->s.endtime - curr->timer->s.starttime;
            //no longer just adds it up, this is for the special purpose full execution timer
                //totalTimes[0] +=curr->timer->s.endtime - curr->timer->s.starttime;
                }
        
        curr = curr->next;
    }
    totalTimes[0] = fullExecTimer.endtime - fullExecTimer.starttime;

}

//assumes simpleNameTally was already called (once) to add up timers
void simpleNamePrint() {
    struct timer_name_tree_node * curr = &root;
    while (curr != NULL) { //still unique names to be checked
        if (checkName(curr->string, rootStr, 0) != 0) {// if the string isn't empty
        printf("Timer [%s]: \t %lu\n", curr->string, curr->times[0]);
        } else {
        printf("Unnamed Timers: \t %lu\n", curr->times[0]);
        }
        printf("\tD2H:    \t %lu\n", curr->times[1]);
        printf("\tH2D:    \t %lu\n", curr->times[2]);
        printf("\tD2D:    \t %lu\n", curr->times[3]);
        printf("\tKernel: \t %lu\n", curr->times[4]);
        printf("\tHost:   \t %lu\n", curr->times[5]);
        printf("\tDual:   \t %lu\n", curr->times[6]);
        curr = curr->next;
    }
}
/*
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
}*/

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
    return (void *) - 1;
}
//only returns a composed timer with events matching both e1 and e2, in either order

void * getDualTimePtr(cl_event e1, cl_event e2) {
    struct timer_group_mem * curr = head.next;
    while (curr != 0) {
        //if composed, will be a negative value
        if (curr->timer->s.type < 0) {
            if ((curr->timer->c.event[0] == e1 \
                    && curr->timer->c.event[1] == e2) \
                    || (curr->timer->c.event[0] == e2 \
                    && curr->timer->c.event[1] == e1)) \
                return (void *) curr->timer;
        }
        curr = curr->next;
    }
    return (void *) - 1;
}

//simply adds timer t to the end of the list

void addTimer(union ocdInternalTimer * t) {
    if (head.next == NULL) { //no members
        tail = &head; //reset tail, just incase
    } else {
        while (tail->next != NULL) { //slide tail to the end
            tail = tail->next;
        }
    }
    struct timer_group_mem * temp_wrap = (struct timer_group_mem *) calloc(sizeof (struct timer_group_mem), 1);

    temp_wrap->next = NULL;
    temp_wrap->timer = t;
    tail->next = temp_wrap;
    tail = temp_wrap;
}

//irreversible! Only call immediately before freeing the timer!

int removeTimer(union ocdInternalTimer * t) {
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


#ifdef TIMER_TEST
//Debug call for checking list construction

void walkList() {
    struct timer_group_mem * curr = head.next;
    fprintf(stderr, "Walking list starting at [%lx]--[%lx]\n", (unsigned long) &head, (unsigned long) head.timer);
    while (curr != 0) {
        fprintf(stderr, "\t[%lx]--[%lx]\n", (unsigned long) curr, (unsigned long) curr->timer->s.event);
        curr = curr->next;
    }
}



//Debug call for checking tree construction recursively

void dfsWalkTree(struct timer_name_tree_node * r, int depth) {
    int i;
    struct timer_name_tree_node * curr = r;
    while (curr != NULL) {
        for (i = depth; i > 0; i--) printf(" ");
        printf("[%s]\n", &curr->string[1]);
        //always descend, the NULL check will kick us back up as needed
        dfsWalkTree(curr->child, depth + curr->len);
        curr = curr->next;
    }
}




//Tests a slew of random, but valid names on the tree or list
//Reuses the tree data type for convenience, but creates a one-way ring queue of generated names
//Names are assigned a random queue longevity to ensure some repetition.
//Assigns random times to each fake event added by abusing the timer type and sidestepping START/END
//USE AT YOUR OWN PERIL DURING A REAL RUN, IT MAY BREAK THE REAL TIMERS
// AND IT *WILL* ADD MASSIVE POSTPROCESSING TIME

void testRandomNames() {
    struct timer_name_tree_node * qhead, *qtail, *temp;
    qhead = qtail = (struct timer_name_tree_node *) calloc(sizeof (struct timer_name_tree_node), 1);
    qhead->next = qhead;
    int i;
    long long ncount = 0, ccount=0;
    srand(time(NULL));
    //generate the names and put them on the queue
    for (i = 0; i < TIMER_TEST_NAME_COUNT; i++) {
        int j = (rand() % TIMER_TEST_MAX_LENG) + 1;
        temp = (struct timer_name_tree_node *) calloc(sizeof (struct timer_name_tree_node), 1);
        temp->string = (char *) calloc(sizeof (char), j + 1);
        temp->len = j; //set the length of the new node
        //fill its name with any standard characters
        for (j--; j >= 0; j--) temp->string[j] = (rand() % 95) + 32;
        temp->tcount = (rand() % TIMER_TEST_MAX_LIFE) + 1; //give it some longevity
        //insert it
        temp->next = qtail->next;
        qtail->next = temp;
        qtail = temp;
    }
    struct timer_name_tree_node * curr = qhead;
    //as long as there's still names on the queue
    while (qhead->next != qhead) {
        if (curr->next == qhead) {
            if (curr->next->next == qhead) continue; //empty!
            curr = curr->next;
        }
        //make a new timer
        struct ocdTimer * timer = (struct ocdTimer *) calloc(sizeof (struct ocdTimer), 1);
        //add it's name data
        timer->type = (enum timer_types)1<<(rand()%5);
        timer->name = curr->next->string;
        timer->nlen = curr->next->len;
        //give it some fake times
        timer->starttime = 0;//rand();
        timer->endtime = 1;// timer->starttime + rand();
        //add it to the global timer list
        addTimer((union ocdInternalTimer *) timer);
        curr->next->tcount--; //reduce the longevity
        if (curr->next->tcount == 0) {
            //remove dead names
            temp = curr->next;
                    curr->next = temp->next;
            free(temp);
        } else { // the node stays alive, move along the ring
            curr = curr->next;
        }
    }
}
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
#define OCD_INIT {\
gettimeofday(&fullExecTimer.timer, NULL);\
fullExecTimer.starttime = 1000 * (fullExecTimer.timer.tv_sec*1000000L + fullExecTimer.timer.tv_usec);\
}

//and absolutely everything needed to finalize them
// performs timer aggregation and printing
// deconstructs timer list and name tree/list
#define OCD_FINISH {\
gettimeofday(&fullExecTimer.timer, NULL);\
fullExecTimer.endtime = 1000 * (fullExecTimer.timer.tv_sec*1000000L + fullExecTimer.timer.tv_usec);\
    simpleNameTally();\
    OCD_PRINT_TIMERS\
    simpleNamePrint();\
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
//assumes must be TIMER_COMPOSED
//#define INI_DUAL_TIMER(a, b) {void * ptr = getDualTimePtr(a, b); \
                         if (ptr == (void *) -1) {\
                                 struct internalComposedTimer * temp = (struct internalComposedTimer*) malloc(sizeof(struct internalComposedTimer)); \
                                 temp->type = TIMER_DUAL;\
                                 temp->starttime = 0;\
                                 temp->endtime = 0;\
                                 temp->event[0] = a;\
                                 temp->event[1] = b;\
                                 addTimer((union internalTimer *)temp);\
                         } else { \
                 		fprintf(stderr, "Timer Error: events [%lx] and [%lx] have already initialized timer [%lx]!\n", (unsigned long) a, (unsigned long) b, (unsigned long) ptr); \
                  	}}

//Destroys a timer permanently - Deprecated
//#define DEST_DUAL_TIMER(t) {void * ptr = getDualTimePtr(a, b); \
                       if(ptr != (void *) -1) { \
                       		removeTimer((union internalTimer *) ptr); \
                       		free(ptr); \
                       }}
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

