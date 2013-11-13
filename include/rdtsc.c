#include "rdtsc.h"

cl_event ocdTempEvent;
#ifdef ENABLE_TIMER
cl_ulong startTime, endTime;

struct ocdTimer * ocdTempTimer;
struct ocdDualTimer * ocdTempDualTimer;
struct ocdHostTimer * ocdTempHostTimer;
struct ocdHostTimer fullExecTimer = {OCD_TIMER_HOST, NULL, 0, 0, 0, {0, 0}};

struct timer_group_mem head = {NULL, NULL, NULL};
struct timer_group_mem * tail;

char rootStr[1] = { (char)0};
cl_ulong rootTimes[7] = {0, 0, 0, 0, 0, 0, 0};
cl_ulong totalTimes[7] = {0, 0, 0, 0, 0, 0, 0};

struct timer_name_tree_node  root = {
    rootStr, 0, NULL, NULL, &head, 0, rootTimes
}; //sentinel

//linear search of the Name List.
//returns a pointer to the correct time array, or -1 if none exists yet
//rather inefficient if many names are used, but the tree will take care of
// speeding lookups, and we'll switch to alpha sort by default as a sideffect
void * checkSimpleNameList(const char * s, int len) {
    struct timer_name_tree_node * curr = root.next;
    while (curr != NULL) { //still unique names to be checked
        if (strcmp(s, curr->string) == 0) {
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
void simpleNameTally()
{
	void * time;
    struct timer_group_mem * curr = head.next;
    while (curr != NULL)
    {
        if (curr->timer->s.name != NULL)
        {
            time = checkSimpleNameList(curr->timer->s.name, curr->timer->s.nlen);
            if (time == (void *)-1)
            {
                //initialize a new name list node
                atail->next = (struct timer_name_tree_node *) calloc(sizeof (struct timer_name_tree_node), 1);
                atail = atail->next;
                atail->next = NULL;
                atail->len = curr->timer->s.nlen;
                atail->string = curr->timer->s.name;
                atail->times = (cl_ulong *) calloc(sizeof(cl_ulong), 7);
                time = (void *)atail->times;
            }
        }
        else
        {
            time = (void *)root.times;
        }
		if (curr->timer->s.endtime > curr->timer->s.starttime)
		{
			switch (curr->timer->s.type)
			{
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
		}
        curr = curr->next;
    }
    totalTimes[0] = fullExecTimer.endtime - fullExecTimer.starttime;
}


//assumes simpleNameTally was already called (once) to add up timers
//now culls off zero-value timers
void simpleNamePrint()
{
    struct timer_name_tree_node * curr = &root;
    while (curr != NULL)
    { //still unique names to be checked
        if (curr->times[0] > 0)
        {
        	if (strcmp(curr->string, rootStr) != 0) // if the string isn't empty
        		printf("Timer [%s]: \t %llu\n", curr->string, curr->times[0]);
        	else
        		printf("Unnamed Timers: \t %llu\n", curr->times[0]);

			if (curr->times[1] > 0) printf("\tD2H:    \t %llu\n", curr->times[1]);
			if (curr->times[2] > 0) printf("\tH2D:    \t %llu\n", curr->times[2]);
			if (curr->times[3] > 0) printf("\tD2D:    \t %llu\n", curr->times[3]);
			if (curr->times[4] > 0) printf("\tKernel: \t %llu\n", curr->times[4]);
			if (curr->times[5] > 0) printf("\tHost:   \t %llu\n", curr->times[5]);
			if (curr->times[6] > 0) printf("\tDual:   \t %llu\n", curr->times[6]);
        }
        curr = curr->next;
    }
}

//Zeroes all timers in name list (e.g., to be used when timing multiple executions without the need to re-allocate/init. We should print before calling.)
void resetNameList()
{
	int ii;
	struct timer_name_tree_node * temp, * curr = root.next;
	while (curr != NULL)
	{
		for(ii=0; ii<7; ii++)
			curr->times[ii] = 0;
		curr=curr->next;
	}
}

//chews up the timer list from head to tail, deallocating all nodes
void destTimerList() {
    struct timer_group_mem * temp, * curr = head.next;
    temp = curr;
    //make sure we can't try to do another cleanup
    head.next = NULL;
    while (curr != NULL) {
        //slide the window
        curr = curr->next;
        //cleanup behind
        if (temp != NULL) {
            if (temp->timer !=NULL) free(temp->timer);
            free(temp);
        }
        //catch up
        temp = curr;
    }
}

//chews up the simpleNameList from root to atail, deallocating all nodes
void destNameList() {
    struct timer_name_tree_node * temp, * curr = root.next;
    temp = curr;
    //make sure we can't try to do another cleanup
    root.next = NULL;
    while (curr != NULL) {
        curr=curr->next;
        if (temp !=NULL) {
            if (temp->times != NULL) free(temp->times);
            free(temp);
        }
        temp = curr;
    }
}

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

#endif //TIMER_TEST

#endif //ENABLE_TIMER