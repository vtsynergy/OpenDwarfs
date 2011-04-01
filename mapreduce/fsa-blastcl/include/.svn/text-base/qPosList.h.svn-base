#ifndef _qPosList_
#define _qPosList_

#ifdef __MY_EXTERN_C__
extern "C" {
#endif

struct queryPosition
{
	uint2 queryPosition;
    struct codeword* codewords;
};

struct codeword
{
	int4 codeword;
	struct codeword* next;
};

extern struct memSingleBlock* qPosList_qPosLists;
extern int4 qPosList_numQPosLists;
extern int4 qPosList_maxQPosLists;

// Initialize for the construction of a new query position list
void qPosList_initialize(int4 maxNumLists);

// Reset between query position list constructions
void qPosList_reset();

// Free the structure
void qPosList_free();

// Print4 the list
void qPosList_print();

// Add a list to the initial collection of query positions lists
void qPosList_addList(int2* queryPositions, int2 numQueryPositions, int4 codeword);

// Process the initial collection of query position lists in order from
// shortest to longest
void qPosList_processLists();

#ifdef __MY_EXTERN_C__
}
#endif

#endif
